import pytest
import hashlib
from eacbp.orchestrator.resources import pin_resource_files
from eacbp.orchestrator.checkpoint import fingerprint


def test_resource_content_changes_resume_signature(tmp_path):
    resource = tmp_path / "network.json"
    resource.write_text("first", encoding="utf-8")
    before = pin_resource_files({"network_path": str(resource), "species": "human"})
    resource.write_text("second", encoding="utf-8")
    after = pin_resource_files({"network_path": str(resource), "species": "human"})
    assert fingerprint(before) != fingerprint(after)
    assert before["network_path"] == after["network_path"]


def test_explicit_missing_resource_is_not_silently_ignored(tmp_path):
    with pytest.raises(FileNotFoundError):
        pin_resource_files({"model_path": str(tmp_path / "absent.pkl")})


def test_output_destination_does_not_change_input_signature(tmp_path):
    output = tmp_path / "corrected.h5"
    before = pin_resource_files({"output_path": str(output)})
    output.write_bytes(b"generated output")
    after = pin_resource_files({"output_path": str(output)})
    assert before == after == {"output_path": str(output.resolve())}


def test_executable_and_checkpoint_contents_change_resume_signature(tmp_path):
    executable = tmp_path / "cellbender"
    checkpoint = tmp_path / "checkpoint.pt"
    executable.write_bytes(b"cellbender-v1")
    checkpoint.write_bytes(b"checkpoint-v1")
    parameters = {
        "executable": str(executable),
        "extra_args": ["--epochs", "10", "--checkpoint", str(checkpoint)],
        "run_cwd": str(tmp_path),
    }

    before = pin_resource_files(parameters)
    assert before["executable"] == str(executable.resolve())
    assert before["external_resource_sha256"] == {
        "executable": hashlib.sha256(b"cellbender-v1").hexdigest(),
        "checkpoint": hashlib.sha256(b"checkpoint-v1").hexdigest(),
    }

    checkpoint.write_bytes(b"checkpoint-v2")
    after_checkpoint = pin_resource_files(parameters)
    assert fingerprint(before) != fingerprint(after_checkpoint)
    assert after_checkpoint["external_resource_sha256"]["executable"] == before["external_resource_sha256"]["executable"]

    executable.write_bytes(b"cellbender-v2")
    after_executable = pin_resource_files(parameters)
    assert fingerprint(after_checkpoint) != fingerprint(after_executable)


def test_checkpoint_argument_must_name_an_existing_local_file(tmp_path):
    executable = tmp_path / "cellbender"
    executable.write_bytes(b"cellbender")
    with pytest.raises(FileNotFoundError):
        pin_resource_files({
            "executable": str(executable),
            "cli_args": ["--checkpoint", str(tmp_path / "missing.pt")],
            "run_cwd": str(tmp_path),
        })
