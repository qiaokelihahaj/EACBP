from pathlib import Path, PurePosixPath, PureWindowsPath

import pandas as pd
import pytest

from scripts.verify_public_kang import (
    MARKER_GENES,
    PROJECT_ROOT,
    _apply_source_integrity,
    _marker_direction_report,
    _parse_lock_packages,
    _sha256_file,
    _validate_external_path,
    _validate_external_work_dir,
)


def test_public_kang_lock_parser_reads_exact_pins_and_ignores_hash_lines():
    packages = _parse_lock_packages(
        "anndata==0.12.19 \\\n"
        "    --hash=sha256:abc\n"
        "pydeseq2==0.5.4 \\\n"
        "    --hash=sha256:def\n"
        "# via anndata\n"
    )

    assert packages == {"anndata": "0.12.19", "pydeseq2": "0.5.4"}


def test_public_kang_marker_check_is_direction_only_and_rejects_wrong_sign():
    table = pd.DataFrame(
        {
            "gene": MARKER_GENES,
            "log2_fold_change": [0.1] * len(MARKER_GENES),
            # Large p/q values remain a sign pass; this is not a significance gate.
            "p_value": [0.9] * len(MARKER_GENES),
            "fdr_q_value": [0.9] * len(MARKER_GENES),
        }
    )
    report = _marker_direction_report(table)
    assert report["status"] == "passed"
    assert report["expected_direction"] == "positive log2 fold change for stim versus ctrl"

    table.loc[table["gene"] == "MX1", "log2_fold_change"] = -0.1
    report = _marker_direction_report(table)
    assert report["status"] == "failed"
    assert next(row for row in report["markers"] if row["gene"] == "MX1")["status"] == "direction_not_positive"


def test_public_kang_marker_check_reports_missing_markers_as_failure():
    table = pd.DataFrame({"gene": ["ISG15"], "log2_fold_change": [0.5]})

    report = _marker_direction_report(table)

    assert report["status"] == "failed"
    assert next(row for row in report["markers"] if row["gene"] == "IFI6")["status"] == "missing_gene"


def test_kang_work_dir_must_be_absolute_and_outside_the_repository(tmp_path):
    with pytest.raises(ValueError, match="absolute path"):
        _validate_external_work_dir(Path("relative-output"))

    repository_output = Path(PROJECT_ROOT) / "outputs"
    with pytest.raises(ValueError, match="outside the repository"):
        _validate_external_work_dir(repository_output)

    external_output = tmp_path / "validation"
    if external_output.is_relative_to(PROJECT_ROOT):
        with pytest.raises(ValueError, match="outside the repository"):
            _validate_external_work_dir(external_output)
    elif isinstance(tmp_path, PureWindowsPath) and tmp_path.drive.casefold() == "c:":
        with pytest.raises(ValueError, match="C: drive"):
            _validate_external_work_dir(external_output)
    else:
        assert _validate_external_work_dir(external_output) == external_output.resolve()

def test_external_path_rules_cover_windows_drives_and_linux_mount_without_io():
    windows_repo = PureWindowsPath("F:/workspace/EACBP")
    with pytest.raises(ValueError, match="absolute path"):
        _validate_external_path(PureWindowsPath("relative-output"), project_root=windows_repo)
    with pytest.raises(ValueError, match="outside the repository"):
        _validate_external_path(windows_repo / "outputs", project_root=windows_repo)
    with pytest.raises(ValueError, match="C: drive"):
        _validate_external_path(PureWindowsPath("C:/validation/kang"), project_root=windows_repo)
    c_drive_repo = PureWindowsPath("C:/workspace/EACBP")
    with pytest.raises(ValueError, match="outside the repository"):
        _validate_external_path(c_drive_repo / "outputs", project_root=c_drive_repo)

    external_windows = PureWindowsPath("F:/validation/kang")
    assert _validate_external_path(external_windows, project_root=PureWindowsPath("C:/workspace/EACBP")) == external_windows

    linux_validation = PurePosixPath("/validation")
    assert _validate_external_path(linux_validation, project_root=PurePosixPath("/workspace/EACBP")) == linux_validation


def test_native_work_dir_validator_uses_current_platform_paths():
    if isinstance(PROJECT_ROOT, PureWindowsPath):
        with pytest.raises(ValueError, match="C: drive"):
            _validate_external_work_dir(Path("C:/validation/kang"))
    else:
        linux_validation = Path("/validation")
        assert _validate_external_work_dir(linux_validation) == linux_validation.resolve()


def test_kang_source_integrity_records_hashes_and_fails_on_mutation(tmp_path):
    source_path = tmp_path / "source.h5ad"
    source_path.write_bytes(b"original source bytes")
    sha256_before = _sha256_file(source_path)
    summary = {"source": {"expected_sha256": sha256_before}, "overall_status": "passed"}

    unchanged = _apply_source_integrity(summary, source_path, sha256_before)

    assert unchanged["unchanged"] is True
    assert summary["source"]["sha256_before"] == sha256_before
    assert summary["source"]["sha256_after"] == sha256_before
    assert summary["overall_status"] == "passed"

    source_path.write_bytes(b"modified source bytes")
    changed_summary = {"source": {}, "overall_status": "passed"}
    changed = _apply_source_integrity(changed_summary, source_path, sha256_before)

    assert changed["unchanged"] is False
    assert changed["sha256_after"] == _sha256_file(source_path)
    assert changed_summary["source"]["sha256_before"] == sha256_before
    assert changed_summary["source"]["sha256_after"] != sha256_before
    assert changed_summary["overall_status"] == "failed"
    assert "changed during the run" in changed_summary["source"]["integrity_failure"]
