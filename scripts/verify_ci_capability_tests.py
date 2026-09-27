"""Run one optional capability suite and fail if any selected test is skipped."""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CAPABILITY_TESTS = {
    "standard": [
        "tests/test_real_methods.py",
        "tests/test_standard_study_script.py",
        "-k",
        "not dpt_to_cellrank_is_connected_in_study_plan",
    ],
    "advanced-statistics": [
        "tests/test_advanced_statistics.py",
        "tests/test_advanced_statistics_contrast_alpha.py",
        "tests/test_advanced_statistics_sparse.py",
        "tests/test_advanced_pipeline.py",
        "tests/test_inference_end_to_end.py",
        "tests/test_inference_audit_contract.py",
    ],
    "advanced-qc": [
        "tests/test_advanced_qc.py",
        "tests/test_celltypist_real.py",
        "tests/test_qc_pipeline.py",
    ],
    "communication": ["tests/test_advanced_communication.py"],
    "fate": [
        "tests/test_cellrank.py",
        "tests/test_real_methods.py::test_dpt_to_cellrank_is_connected_in_study_plan",
    ],
}


def _test_cases(report_path: Path) -> list[ET.Element]:
    try:
        root = ET.parse(report_path).getroot()
    except (ET.ParseError, OSError) as exc:
        raise RuntimeError(f"could not read pytest JUnit report: {exc}") from exc
    return list(root.iter("testcase"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capability", choices=sorted(CAPABILITY_TESTS))
    args = parser.parse_args(argv)

    with tempfile.TemporaryDirectory(prefix="eacbp-capability-tests-") as temp_dir:
        report_path = Path(temp_dir) / "junit.xml"
        command = [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            f"--junitxml={report_path}",
            f"--basetemp={Path(temp_dir) / 'pytest'}",
            *CAPABILITY_TESTS[args.capability],
        ]
        completed = subprocess.run(command, cwd=PROJECT_ROOT, check=False)
        if completed.returncode != 0:
            return completed.returncode

        try:
            cases = _test_cases(report_path)
        except RuntimeError as exc:
            print(f"Capability lane {args.capability}: {exc}", file=sys.stderr)
            return 1

        if not cases:
            print(
                f"Capability lane {args.capability} collected no tests.",
                file=sys.stderr,
            )
            return 1

        skipped = [case for case in cases if case.find("skipped") is not None]
        if skipped:
            print(
                f"Capability lane {args.capability} skipped {len(skipped)} "
                "selected test(s); optional dependencies must be installed and "
                "these tests must run.",
                file=sys.stderr,
            )
            for case in skipped:
                name = case.get("name", "<unnamed>")
                classname = case.get("classname", "<unknown>")
                reason = case.find("skipped").get("message", "no skip reason")
                print(f"  {classname}::{name}: {reason}", file=sys.stderr)
            return 1

        print(
            f"Capability lane {args.capability} ran {len(cases)} test case(s) "
            "with no skips."
        )
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
