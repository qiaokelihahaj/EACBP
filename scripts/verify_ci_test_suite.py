"""Run the explicit core or complete optional-dependency CI test suite."""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from ci_test_suites import PROJECT_ROOT, load_manifest, verify_junit_report


def _run_suite(suite: str, report_path: Path | None) -> int:
    manifest = load_manifest()
    with tempfile.TemporaryDirectory(prefix=f"eacbp-{suite}-tests-") as temp_dir:
        temporary_root = Path(temp_dir)
        output_report = report_path or temporary_root / "junit.xml"
        output_report.parent.mkdir(parents=True, exist_ok=True)
        if suite == "core":
            references = manifest["core"]
            selection_options = [
                f"--deselect={nodeid}" for nodeid in manifest.get("core_deselect", [])
            ]
            allowed_skips = {
                nodeid: reason
                for nodeid, reason in manifest["expected_skips"].items()
                if nodeid.startswith("tests/test_process_supervision.py::")
            }
        else:
            references = ["tests"]
            selection_options = []
            allowed_skips = manifest["expected_skips"]

        command = [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            f"--junitxml={output_report}",
            f"--basetemp={temporary_root / 'pytest'}",
            *selection_options,
            *references,
        ]
        completed = subprocess.run(command, cwd=PROJECT_ROOT, check=False)
        if completed.returncode != 0:
            return completed.returncode

        try:
            count = verify_junit_report(output_report, allowed_skips)
        except RuntimeError as exc:
            print(f"CI {suite} suite: {exc}", file=sys.stderr)
            return 1
        print(f"CI {suite} suite passed: {count} test case(s).")
        return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", choices=("core", "integration"))
    parser.add_argument(
        "--junitxml",
        type=Path,
        help="also save the pytest JUnit report at this path",
    )
    args = parser.parse_args(argv)
    report_path = args.junitxml.resolve() if args.junitxml else None
    return _run_suite(args.suite, report_path)


if __name__ == "__main__":
    raise SystemExit(main())
