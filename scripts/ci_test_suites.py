"""Shared manifest and JUnit checks for the CI test lanes."""

from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = PROJECT_ROOT / ".github" / "ci-test-suites.json"


def load_manifest() -> dict[str, Any]:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    validate_manifest(manifest)
    return manifest


def _module_path(reference: str) -> str:
    return reference.partition("::")[0]


def validate_manifest(manifest: dict[str, Any]) -> None:
    """Fail closed when a test module is not assigned to a CI lane."""
    core = manifest.get("core")
    core_deselect = manifest.get("core_deselect", [])
    capabilities = manifest.get("capabilities")
    integration_only = manifest.get("integration_only")
    if not isinstance(core, list) or not isinstance(capabilities, dict):
        raise RuntimeError("CI test manifest must define core and capabilities")
    if not isinstance(core_deselect, list):
        raise RuntimeError("CI test manifest core_deselect must be a list")
    if not isinstance(integration_only, list):
        raise RuntimeError("CI test manifest must define integration_only")

    references = [*core, *core_deselect, *integration_only]
    for suite in capabilities.values():
        if not isinstance(suite, dict) or not isinstance(suite.get("tests"), list):
            raise RuntimeError("each capability must define a tests list")
        if not isinstance(suite.get("deselect"), list):
            raise RuntimeError("each capability must define a deselect list")
        references.extend(suite["tests"])
        references.extend(suite["deselect"])

    declared_modules = {_module_path(reference) for reference in references}
    available_modules = {
        path.relative_to(PROJECT_ROOT).as_posix()
        for path in (PROJECT_ROOT / "tests").glob("test_*.py")
    }
    missing = sorted(available_modules - declared_modules)
    unknown = sorted(declared_modules - available_modules)
    if missing or unknown:
        details = []
        if missing:
            details.append(f"unassigned test modules: {', '.join(missing)}")
        if unknown:
            details.append(f"missing test modules: {', '.join(unknown)}")
        raise RuntimeError("; ".join(details))

    if len(core) != len(set(core)):
        raise RuntimeError("core test module list contains duplicates")
    if len(core_deselect) != len(set(core_deselect)):
        raise RuntimeError("core deselection list contains duplicates")
    for nodeid in core_deselect:
        if "::" not in nodeid or _module_path(nodeid) not in set(core):
            raise RuntimeError(f"core deselection must target a selected module: {nodeid}")
    if set(core) & set(integration_only):
        raise RuntimeError("core and integration_only test modules overlap")

    capability_references = {
        reference
        for suite in capabilities.values()
        for reference in suite["tests"]
    }
    for suite in capabilities.values():
        for excluded in suite["deselect"]:
            if excluded not in capability_references:
                raise RuntimeError(
                    f"deselected test is not selected by a capability lane: {excluded}"
                )

    expected_skips = manifest.get("expected_skips")
    if not isinstance(expected_skips, dict):
        raise RuntimeError("CI test manifest must define expected_skips")
    for nodeid, reason in expected_skips.items():
        module = _module_path(nodeid)
        if module not in available_modules or "::" not in nodeid or not reason:
            raise RuntimeError(f"invalid expected skip entry: {nodeid}")


def selected_module_path(reference: str) -> Path:
    return PROJECT_ROOT / _module_path(reference)


def verify_junit_report(report_path: Path, expected_skips: dict[str, str]) -> int:
    try:
        cases = list(ET.parse(report_path).getroot().iter("testcase"))
    except (ET.ParseError, OSError) as exc:
        raise RuntimeError(f"could not read pytest JUnit report: {exc}") from exc
    if not cases:
        raise RuntimeError("pytest JUnit report contains no test cases")

    unexpected: list[tuple[str, str]] = []
    observed_expected: set[str] = set()
    for case in cases:
        skipped = case.find("skipped")
        if skipped is None:
            continue
        module_name = (case.get("classname") or "").rsplit(".", 1)[-1]
        test_name = (case.get("name") or "").split("[", 1)[0]
        nodeid = f"tests/{module_name}.py::{test_name}"
        reason = skipped.get("message", "")
        expected_reason = expected_skips.get(nodeid)
        if expected_reason is None or not reason.endswith(expected_reason):
            unexpected.append((nodeid, reason or "no skip reason"))
        else:
            observed_expected.add(nodeid)

    if unexpected:
        lines = [f"  {nodeid}: {reason}" for nodeid, reason in unexpected]
        raise RuntimeError(
            f"pytest skipped {len(unexpected)} unexpected selected test(s):\n"
            + "\n".join(lines)
        )
    if observed_expected:
        labels = ", ".join(sorted(observed_expected))
        print(f"Allowed documented environment skip(s): {labels}")
    print(f"JUnit report contains {len(cases)} test case(s) with no unexpected skips.")
    return len(cases)
