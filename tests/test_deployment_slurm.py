import json
from pathlib import Path
import subprocess

import pytest

import eacbp.deployment.slurm as deployment_slurm
from eacbp.deployment.slurm import (
    DeploymentError,
    build_submit_command,
    default_batch_script,
    main,
    validate_job_environment,
)


def _job_environment(tmp_path: Path, *, action: str = "run", run_id: str = "study-001") -> dict[str, str]:
    image = tmp_path / "eacbp.sif"
    image.write_bytes(b"test image placeholder")
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    config = tmp_path / "config"
    config.mkdir()
    data = tmp_path / "data"
    data.mkdir()
    (config / "manifest.json").write_text("{}", encoding="utf-8")
    (data / "study.h5ad").write_bytes(b"fixture")
    values = {
        "EACBP_ACTION": action,
        "EACBP_RUN_ID": run_id,
        "EACBP_IMAGE": str(image),
        "EACBP_OUTPUTS": str(outputs),
    }
    if action == "run":
        values.update({
            "EACBP_CONFIG_DIR": str(config),
            "EACBP_MANIFEST": "manifest.json",
            "EACBP_DATA_DIR": str(data),
            "EACBP_DATA": "study.h5ad",
        })
    else:
        (outputs / run_id).mkdir()
    return values


def test_build_submit_command_uses_argument_list_and_resource_overrides(tmp_path: Path) -> None:
    command = build_submit_command(
        script=tmp_path / "run.sbatch",
        logs_dir=tmp_path / "study logs",
        cwd=tmp_path,
        cpus_per_task=8,
        memory="64G",
        time_limit="1-12:00:00",
        partition="compute",
        account="bio_lab",
    )

    assert command[0] == "sbatch"
    assert "--cpus-per-task=8" in command
    assert "--mem=64G" in command
    assert "--time=1-12:00:00" in command
    assert f"--output={tmp_path / 'study logs'}/%x-%j.out" in command
    assert command[-1] == str(tmp_path / "run.sbatch")
    assert all(";" not in arg for arg in command)


def test_default_batch_script_resolves_to_canonical_template() -> None:
    script = default_batch_script()

    assert script.is_file()
    assert script.name == "run_eacbp_container.sbatch"


def test_build_submit_command_rejects_shell_metacharacter_in_partition(tmp_path: Path) -> None:
    with pytest.raises(DeploymentError, match="partition"):
        build_submit_command(
            script=tmp_path / "run.sbatch",
            logs_dir=tmp_path,
            cwd=tmp_path,
            partition="compute; touch /tmp/unsafe",
        )


def test_validate_run_environment_resolves_inputs_inside_read_only_roots(tmp_path: Path) -> None:
    normalized = validate_job_environment(_job_environment(tmp_path))

    assert normalized["EACBP_DATA"] == "study.h5ad"
    assert normalized["EACBP_MANIFEST"] == "manifest.json"
    assert normalized["EACBP_ACTION"] == "run"


def test_validate_run_environment_rejects_input_path_escape(tmp_path: Path) -> None:
    values = _job_environment(tmp_path)
    values["EACBP_DATA"] = "../outside.h5ad"

    with pytest.raises(DeploymentError, match="relative path contained"):
        validate_job_environment(values)


def test_dry_run_validates_without_creating_logs_or_contacting_slurm(tmp_path: Path, monkeypatch, capsys) -> None:
    values = _job_environment(tmp_path)
    for key, value in values.items():
        monkeypatch.setenv(key, value)
    script = tmp_path / "run.sbatch"
    script.write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    logs = tmp_path / "logs"

    status = main([
        "submit", "--script", str(script), "--logs-dir", str(logs),
        "--cwd", str(tmp_path), "--dry-run",
    ])

    output = json.loads(capsys.readouterr().out)
    assert status == 0
    assert output["scheduler_contacted"] is False
    assert output["action"] == "run"
    assert not logs.exists()


def test_submit_forwards_normalized_environment_to_sbatch(tmp_path: Path, monkeypatch, capsys) -> None:
    values = _job_environment(tmp_path)
    alias_dir = tmp_path / "alias"
    alias_dir.mkdir()
    values["EACBP_IMAGE"] = str(alias_dir / ".." / "eacbp.sif")
    values["EACBP_OUTPUTS"] = str(alias_dir / ".." / "outputs")
    for key, value in values.items():
        monkeypatch.setenv(key, value)
    script = tmp_path / "run.sbatch"
    script.write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    logs = tmp_path / "logs"
    captured: dict[str, object] = {}

    def fake_run(command, *, env=None):
        captured["env"] = env
        return subprocess.CompletedProcess(command, 0, "12345;cluster\n", "")

    monkeypatch.setattr(deployment_slurm, "_require_command", lambda name: name)
    monkeypatch.setattr(deployment_slurm, "_run", fake_run)
    status = main([
        "submit", "--script", str(script), "--logs-dir", str(logs),
        "--cwd", str(tmp_path),
    ])

    output = json.loads(capsys.readouterr().out)
    forwarded = captured["env"]
    assert status == 0
    assert output["job_id"] == "12345"
    assert isinstance(forwarded, dict)
    assert forwarded["EACBP_IMAGE"] == str(tmp_path / "eacbp.sif")
    assert forwarded["EACBP_OUTPUTS"] == str(tmp_path / "outputs")


def test_run_environment_rejects_existing_run_dir(tmp_path: Path) -> None:
    values = _job_environment(tmp_path)
    (Path(values["EACBP_OUTPUTS"]) / values["EACBP_RUN_ID"]).mkdir()

    with pytest.raises(DeploymentError, match="already exists"):
        validate_job_environment(values)


def test_resume_can_rebind_persisted_reference_root(tmp_path: Path) -> None:
    values = _job_environment(tmp_path, action="resume")
    refs = tmp_path / "refs"
    refs.mkdir()
    values["EACBP_REFS_DIR"] = str(refs)

    normalized = validate_job_environment(values)

    assert normalized["EACBP_REFS_DIR"] == str(refs.resolve())
