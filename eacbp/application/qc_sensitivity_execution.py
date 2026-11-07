"""Execute frozen QC branches and publish a detached, verified comparison.

Each run retains its native registry, audit receipts, snapshot and report. The
comparison is a derived analysis and cannot replace any run's scientific status.
"""
from __future__ import annotations

from pathlib import Path
import time
from typing import Any

from eacbp._atomic_json import atomic_write_json
from eacbp.schemas.artifact import ArtifactType
from eacbp.application.qc_sensitivity import (
    compare_qc_sensitivity_runs, qc_branch_config, qc_config_sha256,
)


def collect_qc_run(run_dir: Path, *, source_sha256: str, config, parameters) -> dict[str, Any]:
    from eacbp.application import study_service
    from eacbp.evidence.snapshot import load_study_snapshot

    registry = study_service._make_registry(run_dir / "artifacts")
    snapshot = load_study_snapshot(run_dir / "snapshot.json", artifact_registry=registry)
    actual_config = snapshot.config
    if actual_config.get("autonomous_analysis", {}).get("source_sha256") != source_sha256:
        raise ValueError("Native QC run snapshot has a different raw source")
    if qc_config_sha256(actual_config) != qc_config_sha256(config):
        raise ValueError("Native QC run configuration changed from the frozen scope")
    receipt = {
        "source_sha256": source_sha256, "config_sha256": qc_config_sha256(actual_config),
        "parameters": dict(parameters), "qc_admitted": False,
        "run_status": snapshot.summary.get("status"), "deg_tables": [], "summaries": {},
        "provenance": {"run_dir": str(run_dir), "snapshot": str(snapshot.path),
                       "snapshot_sha256": study_service._hash_file(snapshot.path),
                       "summary_sha256": study_service._hash_file(run_dir / "summary.json"),
                       "config_sha256": study_service._hash_file(run_dir / "config.json"),
                       "report": str(run_dir / "report.md"), "artifacts": []},
    }
    admitted = {context.task_id for context in snapshot.audit_contexts.values()
                if context.overall_passed and context.receipt_present}
    for task in snapshot.task_results:
        if task.task_id not in admitted or task.status.value != "success":
            continue
        if task.capability not in {"qc", "qc_review", "clustering", "deg", "donor_sensitivity", "cluster_markers", "rna_velocity"}:
            continue
        receipt["summaries"][task.task_id] = dict(task.metrics)
        for uri in task.output_artifacts:
            meta = registry.get_metadata(uri)
            receipt["provenance"]["artifacts"].append({
                "uri": uri, "sha256": meta.sha256_hash, "task_id": task.task_id,
                "audit_status": "admitted",
            })
            if task.capability == "qc" and meta.type == ArtifactType.ANNDATA:
                _, payload = registry.get_scdata(uri)
                actual = {key: meta.parameters.get(key) for key in ("min_genes", "max_mito_pct")}
                if actual != parameters:
                    raise ValueError("Native QC artifact parameters differ from the frozen thresholds")
                receipt["parameters"] = actual
                receipt["cell_disposition"] = payload.uns["qc"]["cell_disposition"]
                receipt["qc_admitted"] = True
            elif task.capability == "clustering" and meta.type == ArtifactType.ANNDATA:
                _, payload = registry.get_scdata(uri)
                if "cluster" in payload.obs:
                    receipt["clusters"] = payload.obs["cluster"].copy()
            elif task.capability == "deg" and meta.type == ArtifactType.TABLE:
                _, table = registry.get(uri)
                receipt["deg_tables"].append({"target_cell_type": task.metrics.get("target_cell_type"),
                                              "table": table, "admitted": True})
    return receipt


def _markdown(result) -> str:
    lines = ["# QC sensitivity comparison", "", f"Status: `{result['status']}`; frozen plan: `{result['plan_id']}`.",
             "", result["interpretation"], "", "The primary native analysis retains its original thresholds and status.",
             "", f"Completion scope: {result.get('completion_scope', 'Frozen QC execution and paired comparisons only; scientific robustness is assessed separately')}.",
             "", "| Branch | QC minimum genes / maximum mitochondrial % | Native status | Comparison | Retained cells | Added / removed | Shared-cell ARI |",
             "| :--- | :--- | :--- | :--- | ---: | :--- | ---: |"]
    for item in result["branches"]:
        params = item["parameters"]
        cluster = item.get("clustering", {})
        ari = cluster.get("shared_cell_ari", "unavailable")
        lines.append(f"| {item['branch_id']} | {params['min_genes']} / {params['max_mito_pct']:.6g} | {item.get('run_status', 'unknown')} | {item['status']} | {item.get('branch_retained', 'unavailable')} | {item.get('n_added', 'unavailable')} / {item.get('n_removed', 'unavailable')} | {ari} |")
    reasons = [f"- {item['branch_id']}: {item['reason']}" for item in result["branches"] if item.get("reason")]
    reasons.extend(f"- {item['branch_id']}: {item['status']}: {item['reason']}" for item in result.get("failures", []))
    if reasons:
        lines.extend(["", *reasons])
    lines.extend(["", "## Native run provenance", ""])
    for branch_id, record in result.get("run_provenance", {}).items():
        lines.append(f"- {branch_id}: [native report]({record['report']}); snapshot SHA-256 `{record['snapshot_sha256']}`.")
    for item in result["branches"]:
        if item.get("observed_group_retention"):
            lines.extend(["", f"## {item['branch_id']}: observed population retention", "",
                          "| Observed population | Original | Primary retained | Branch retained |", "| :--- | ---: | ---: | ---: |"])
            for group in item["observed_group_retention"]:
                lines.append(f"| {group['observed_group']} | {group['original_cells']} | {group['primary_retained']} | {group['branch_retained']} |")
        if "deg" in item:
            import json
            lines.extend(["", f"## {item['branch_id']}: differential result sensitivity", "", "```json",
                          json.dumps({key: value for key, value in item["deg"].items() if key != "rows"}, ensure_ascii=False, indent=2), "```"])
    lines.extend(["", "Full gene comparisons, cell dispositions, artifact hashes and unavailable-result reasons are retained in `comparison.json`.",
                  "Shared-cell agreement and repeated differential estimates are dependent sensitivity evidence. They do not establish independent biological replication.",
                  "The derived comparison and branch registries are separate from the primary study export; preserve this entire QC sensitivity directory together with the primary run.", ""])
    return "\n".join(lines)


def execute_qc_sensitivity(*, plan, manifest, config, data: Path, run_dir: Path,
                           baseline_summary, observed_groups=None, started_at=None):
    from eacbp.application import study_service

    # Test doubles and failed runs may have no persisted native snapshot. They
    # supply no admissible sensitivity evidence and must not trigger replays.
    if not (run_dir / "snapshot.json").is_file():
        return baseline_summary
    started_at = time.monotonic() if started_at is None else started_at
    baseline_files = {name: study_service._hash_file(run_dir / name)
                      for name in ("snapshot.json", "summary.json", "report.md", "config.json", "artifacts/.artifact_registry.json")
                      if (run_dir / name).is_file()}
    root = run_dir / "qc_sensitivity"
    root.mkdir()
    atomic_write_json(root / "plan.json", plan)
    plan_file_sha256 = study_service._hash_file(root / "plan.json")
    try:
        baseline = collect_qc_run(run_dir, source_sha256=plan["source_sha256"], config=config,
                                  parameters=plan["primary"]["parameters"])
        branches, failures = {}, []
        for spec in plan["branches"]:
            if spec["planned_status"] != "planned":
                continue
            branch_id = spec["branch_id"]
            if not baseline["qc_admitted"]:
                failures.append({"branch_id": branch_id, "status": "not_run", "reason": "Primary QC was not independently admitted"})
                continue
            if study_service._hash_file(data) != plan["source_sha256"]:
                raise ValueError("Raw source changed after QC planning")
            branch_config = qc_branch_config(plan, branch_id, config)
            branch_manifest = manifest.model_copy(deep=True)
            branch_manifest.study_id += "_qc_" + branch_id
            remaining = manifest.constraints.max_runtime_hours - (time.monotonic() - started_at) / 3600
            if remaining <= 0:
                failures.append({"branch_id": branch_id, "status": "not_run", "reason": "Total analysis runtime budget exhausted"})
                continue
            branch_manifest.constraints.max_runtime_hours = remaining
            branch_dir = root / "runs" / branch_id
            try:
                study_service.run_study(manifest=branch_manifest, config=branch_config, data=data, run_dir=branch_dir)
                branches[branch_id] = collect_qc_run(branch_dir, source_sha256=plan["source_sha256"],
                                                    config=branch_config, parameters=spec["parameters"])
            except Exception as exc:
                failures.append({"branch_id": branch_id, "status": "failed", "reason": f"{type(exc).__name__}: {exc}", "run_dir": str(branch_dir)})
        result = compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches, observed_groups=observed_groups)
        result["failures"] = failures
        result["run_provenance"] = {"primary": baseline["provenance"], **{key: value["provenance"] for key, value in branches.items()}}
        for name, digest in baseline_files.items():
            if study_service._hash_file(run_dir / name) != digest:
                raise ValueError("Primary native evidence changed during QC sensitivity execution")
        # Recheck every registry and payload after the other branches finish;
        # an earlier verified receipt cannot bless a later changed artifact.
        from eacbp.evidence.snapshot import load_study_snapshot
        for record in result["run_provenance"].values():
            path = Path(record["run_dir"]) if record.get("run_dir") else None
            if path is not None:
                load_study_snapshot(path / "snapshot.json", artifact_registry=study_service._make_registry(path / "artifacts"))
                for name in ("snapshot", "summary", "config"):
                    if study_service._hash_file(path / (name + ".json")) != record[name + "_sha256"]:
                        raise ValueError("Native run byte hashes changed after sensitivity collection")
        if study_service._hash_file(root / "plan.json") != plan_file_sha256:
            raise ValueError("Frozen QC plan file changed during execution")
        result["plan_file_sha256"] = plan_file_sha256
        atomic_write_json(root / "comparison.json", result)
        (root / "report.md").write_text(_markdown(result), encoding="utf-8")
        return {**baseline_summary, "qc_sensitivity": {"status": result["status"], "plan_id": plan["plan_id"],
                 "report": str(root / "report.md"), "comparison": str(root / "comparison.json")}}
    except Exception as exc:
        failure = {"status": "failed", "plan_id": plan["plan_id"], "reason": f"{type(exc).__name__}: {exc}",
                   "primary_evidence_preserved": all(study_service._hash_file(run_dir / name) == digest for name, digest in baseline_files.items())}
        atomic_write_json(root / "comparison.json", failure)
        return {**baseline_summary, "qc_sensitivity": {**failure, "comparison": str(root / "comparison.json")}}
