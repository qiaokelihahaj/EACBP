"""Replay the existing Endometrium_MU small-table exports through EACBP.

No SSH, FASTQ downloads, h5ad computation or changes to the source project.
Expected source hashes come from its historical validation index, not from
accepting whatever bytes happen to be at the source paths today.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PureWindowsPath
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eacbp.application.descriptive_review import review_descriptive_bundle, publish_descriptive_review


def verify(project: Path, run_dir: Path) -> dict:
    source_dir = project / "results" / "epithelial_molecular_report_v1"
    validation_bytes = (source_dir / "validation.json").read_bytes()
    validation = json.loads(validation_bytes)
    roles = {
        "cells": "driver_gene_raw_UMI_per_cell.csv",
        "gene_summaries": "gene_raw_UMI_by_mask_sample.csv",
        "module_summaries": "module_deltas_by_group_mask.csv",
        "scoring": "run_manifest.json",
        "sensitivity": "sensitivity_manifest.json",
    }
    payloads, hashes = {}, {}
    for role, name in roles.items():
        matches = [item["sha256"] for item in validation["source_SHA_verified"]
                   if PureWindowsPath(item["path"]).name == name]
        if len(matches) != 1:
            raise ValueError(f"Historical index does not uniquely pin {name}")
        data = (source_dir / "data" / name).read_bytes()
        if hashlib.sha256(data).hexdigest() != matches[0]:
            raise ValueError(f"Historical source SHA256 mismatch: {name}")
        payloads[role], hashes[role] = data, matches[0]
    roles["modules"] = "epithelial_modules_v1_1.json"
    payloads["modules"] = (project / "metadata" / roles["modules"]).read_bytes()
    hashes["modules"] = json.loads(payloads["scoring"])["module_file_sha256"]
    if hashlib.sha256(payloads["modules"]).hexdigest() != hashes["modules"]:
        raise ValueError("Module definitions differ from the historical scoring run")
    # Validate all bytes before creating any review output. Existing runs stay intact.
    run_dir = run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / "sources").mkdir()
    for role, name in roles.items():
        (run_dir / "sources" / name).write_bytes(payloads[role])
    (run_dir / "historical_validation.json").write_bytes(validation_bytes)
    spec = {
        "schema_version": 1, "format": "frozen_scoring_tables_v1", "study_id": "Endometrium_MU",
        "population": "epithelial_unresolved", "reference_sample": "P12K8_con", "query_sample": "P12K8_cKO",
        "min_cells": 20, "design": {"observation_unit": "library", "biological_replication": "unverified"},
        "sources": {role: {"path": "sources/" + name, "sha256": hashes[role]} for role, name in roles.items()},
    }
    bundle = run_dir / "bundle.json"
    bundle.write_text(json.dumps(spec, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    digest = hashlib.sha256(bundle.read_bytes()).hexdigest()
    result = review_descriptive_bundle(bundle, expected_bundle_sha256=digest)
    publish_descriptive_review(result, output=run_dir / "review.json", report=run_dir / "report.md")
    actual = result["checks"]["checks"]
    expected = {"n_cells": 515, "n_raw_genes": 16, "n_masks": 7,
                "n_distinct_mask_cell_sets": 5, "n_module_ids": 19, "n_module_rows_checked": 931}
    mismatches = {key: {"actual": actual.get(key), "expected": value}
                  for key, value in expected.items() if actual.get(key) != value}
    receipt = {"status": "passed" if not mismatches else "failed", "scope": "real_export_replay",
               "bundle_sha256": digest, "historical_validation_sha256": hashlib.sha256(validation_bytes).hexdigest(),
               "source_project": str(project.resolve()), "expected": expected, "observed": actual,
               "mismatches": mismatches, "report": str(run_dir / "report.md"),
               "review_sha256": hashlib.sha256((run_dir / "review.json").read_bytes()).hexdigest(),
               "report_sha256": hashlib.sha256((run_dir / "report.md").read_bytes()).hexdigest()}
    (run_dir / "acceptance.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True, help="new directory; never overwrites an existing run")
    args = parser.parse_args()
    receipt = verify(args.project, args.run_dir)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0 if receipt["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
