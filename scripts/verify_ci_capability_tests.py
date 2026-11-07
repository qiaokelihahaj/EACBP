"""Run one optional capability suite and fail if any selected test is skipped."""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from ci_test_suites import (
    PROJECT_ROOT,
    load_manifest,
    verify_junit_report,
)


def main(argv: list[str] | None = None) -> int:
    manifest = load_manifest()
    capability_tests = manifest["capabilities"]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capability", choices=sorted(capability_tests))
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
            *capability_tests[args.capability]["tests"],
            *[
                f"--deselect={nodeid}"
                for nodeid in capability_tests[args.capability]["deselect"]
            ],
        ]
        completed = subprocess.run(command, cwd=PROJECT_ROOT, check=False)
        if completed.returncode != 0:
            return completed.returncode

        try:
            count = verify_junit_report(report_path, {})
        except RuntimeError as exc:
            print(f"Capability lane {args.capability}: {exc}", file=sys.stderr)
            return 1

        if count == 0:
            return 1

        print(
            f"Capability lane {args.capability} ran {count} test case(s) "
            "with no skips."
        )
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
