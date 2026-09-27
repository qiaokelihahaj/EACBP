"""Run acceptance with ``python -m eacbp.acceptance MANIFEST.json``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from pydantic import ValidationError

from eacbp.acceptance.manifest import load_manifest
from eacbp.acceptance.runner import acceptance_exit_code, run_acceptance, write_result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Validate donor-level pseudobulk inputs and compare EACBP with direct PyDESeq2.")
    parser.add_argument("manifest", help="Path to a typed research acceptance JSON manifest")
    parser.add_argument("--output", help="Write machine-readable JSON result to this path; defaults to stdout")
    args = parser.parse_args(argv)
    try:
        manifest, manifest_dir = load_manifest(args.manifest)
        result = run_acceptance(manifest, manifest_dir)
    except (OSError, json.JSONDecodeError, ValidationError, ValueError) as exc:
        result = {
            "schema_version": 1,
            "statuses": {"engineering": "failed", "method": "not_run"},
            "research_acceptance": {"status": "blocked", "biological_validity_claimed": False},
            "failure": {"type": type(exc).__name__, "message": str(exc)},
        }
    write_result(result, args.output)
    return acceptance_exit_code(result)


if __name__ == "__main__":
    sys.exit(main())
