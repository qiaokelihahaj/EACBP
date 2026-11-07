from __future__ import annotations

import pytest

from scripts.ci_test_suites import load_manifest, validate_manifest, verify_junit_report


def test_manifest_rejects_an_unassigned_test_module():
    manifest = load_manifest()
    manifest["core"].remove("tests/test_adapters.py")

    with pytest.raises(RuntimeError, match="unassigned test modules: tests/test_adapters.py"):
        validate_manifest(manifest)


def test_junit_gate_rejects_unexpected_skip(tmp_path):
    report = tmp_path / "junit.xml"
    report.write_text(
        '<testsuite tests="1"><testcase classname="tests.test_sample" '
        'name="test_case"><skipped message="missing dependency" /></testcase></testsuite>',
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="tests/test_sample.py::test_case: missing dependency"):
        verify_junit_report(report, {})


def test_junit_gate_accepts_an_exact_documented_skip(tmp_path):
    nodeid = "tests/test_process_supervision.py::test_failed_job_assignment_kills_suspended_parent"
    reason = "Windows must assign the job before resuming the process"
    report = tmp_path / "junit.xml"
    report.write_text(
        '<testsuite tests="1"><testcase classname="tests.test_process_supervision" '
        'name="test_failed_job_assignment_kills_suspended_parent">'
        f'<skipped message="{reason}" /></testcase></testsuite>',
        encoding="utf-8",
    )

    assert verify_junit_report(report, {nodeid: reason}) == 1
