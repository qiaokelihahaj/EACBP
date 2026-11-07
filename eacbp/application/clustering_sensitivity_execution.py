"""Prepare frozen raw-count replays and collect their native audited evidence.

The comparison is detached from every study. There is no winner selection,
QC rescue, automatic invocation from analyze, or external-label input.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from eacbp._atomic_json import atomic_write_json
from eacbp.application.clustering_sensitivity import (
    clustering_branch_config, clustering_config_sha256,
    clustering_scientific_config_sha256, compare_clustering_sensitivity_runs,
    feature_axis_fingerprint, freeze_clustering_sensitivity_plan,
    raw_count_fingerprint, validate_clustering_branch_config,
    validate_clustering_sensitivity_plan,
)
from eacbp.orchestrator.checkpoint import execution_environment
from eacbp.schemas.audit import model_to_jsonable, stable_fingerprint
from eacbp.schemas.runtime import RunConfig
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import TaskContract


def _write(path, payload):
    return atomic_write_json(path, payload, dump_kwargs={"indent": 2, "sort_keys": True,
                                                       "ensure_ascii": False, "allow_nan": False}, create_parent=True)


def _read(path):
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("Frozen clustering document must be a JSON object")
    return value


def _contained(path, root):
    value = Path(path).expanduser().resolve()
    value.relative_to(root)
    return value


def _spec(plan, branch_id):
    return next(item for item in [plan["primary"], *plan["branches"]] if item["branch_id"] == branch_id)


def _counts(data):
    if "counts" not in data.layers:
        raise ValueError("Native prepared and replay inputs must preserve layers.counts")
    return data.layers["counts"]


def _role_values(frame, columns):
    return {column: [None if pd.isna(value) else str(value).strip() for value in frame[column]]
            for column in columns}


def _verify_subset(data, spec):
    if data.obs.index.astype(str).tolist() != spec["selected_ids"]:
        raise ValueError("Derived raw input cell membership/order differs from frozen mask")
    if raw_count_fingerprint(_counts(data)) != spec["raw_counts_sha256"]:
        raise ValueError("Derived raw input counts/dtype differ from the parent source subset")
    if feature_axis_fingerprint(data.var, "var.gene_name") != spec["feature_axis_sha256"]:
        raise ValueError("Derived raw input changed the complete feature universe")


def prepare_clustering_sensitivity(*, manifest, config, data, output_dir) -> dict[str, Any]:
    """Write all derived inputs and frozen configurations before any native run.

    Callers choose the storage directory. A separate execution call permits a
    cross-dataset protocol to pin both preparations before either is executed.
    """
    from eacbp.application import study_service
    from eacbp.application.autonomous import observe_input, prepare_input

    root = Path(output_dir).expanduser().resolve()
    if root.exists():
        raise ValueError("Clustering sensitivity output directory already exists")
    source = Path(data).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError("Clustering sensitivity requires an existing source H5AD")
    baseline = RunConfig.from_mapping(study_service.normalise_config_paths(config, base_dir=Path.cwd())).to_mapping()
    source_sha256 = study_service._hash_file(source)
    payload = study_service._import_h5ad(source)
    settings = baseline.get("autonomous_analysis", {})
    observation = observe_input(payload, question=settings.get("question"))
    if study_service._hash_file(source) != source_sha256:
        raise ValueError("Raw source changed during pre-execution observation")
    plan = freeze_clustering_sensitivity_plan(payload=payload, observations=observation,
                                              source_sha256=source_sha256, baseline_config=baseline)
    canonical, _ = prepare_input(payload, settings)
    del payload
    root.parent.mkdir(parents=True, exist_ok=True)
    root.mkdir(exist_ok=False)
    _write(root / "plan.json", plan)
    environment = execution_environment()
    records = []
    for spec in [plan["primary"], *plan["branches"]]:
        branch_id = spec["branch_id"]
        if spec["planned_status"] != "planned":
            records.append({"branch_id": branch_id, "planned_status": spec["planned_status"], "reason": spec["reason"]})
            continue
        branch_source = source
        if spec["source_kind"] == "derived_subset":
            branch_source = root / "inputs" / (branch_id + ".h5ad")
            branch_source.parent.mkdir(exist_ok=True)
            selected = canonical.subset_obs(np.isin(np.arange(canonical.n_obs), spec["selected_positions"]))
            _verify_subset(selected, spec)
            selected.to_anndata().write_h5ad(branch_source)
            del selected
            restored = study_service._import_h5ad(branch_source)
            _verify_subset(restored, spec)
            if observe_input(restored, question=settings["question"]) != spec["expected_observations"]:
                raise ValueError("Derived H5AD roundtrip changed its frozen observations")
            del restored
        branch_source_sha256 = study_service._hash_file(branch_source)
        branch_config = clustering_branch_config(plan, branch_id, baseline, source_sha256=branch_source_sha256)
        branch_manifest = manifest.model_copy(deep=True)
        identity = hashlib.sha256((str(manifest.study_id) + plan["plan_sha256"] + branch_id).encode()).hexdigest()[:20]
        branch_manifest.study_id = "clustering_sensitivity_" + identity
        branch_manifest.data.raw_artifact_uri = None
        branch_manifest.data.has_spatial_coordinates = bool(spec["expected_observations"]["spatial"]["available"])
        branch_observation = spec["expected_observations"]
        donor_column = branch_observation["obs_column_map"].get("donor")
        branch_manifest.experimental_design.biological_unit = "donor" if donor_column else "unknown"
        branch_manifest.experimental_design.total_samples = (int(canonical.obs.iloc[spec["selected_positions"]][donor_column].nunique())
                                                              if donor_column else 0)
        branch_manifest.experimental_design.batches = list(branch_observation["batches"])
        branch_manifest.biological_design.conditions = list(branch_observation["conditions"])
        role_mappings = [observation["obs_column_map"], branch_observation["obs_column_map"]]
        role_columns = sorted({column for mapping in role_mappings for role, column in mapping.items()
                               if column is not None and role != "cell_type" and column in canonical.obs}
                              | {role for mapping in role_mappings for role, column in mapping.items()
                                 if column is not None and role != "cell_type" and role in canonical.obs})
        role_values = _role_values(canonical.obs.iloc[spec["selected_positions"]], role_columns)
        config_path = root / "frozen_configs" / (branch_id + ".json")
        manifest_path = root / "frozen_manifests" / (branch_id + ".json")
        _write(config_path, branch_config)
        _write(manifest_path, branch_manifest.model_dump(mode="json"))
        records.append({"branch_id": branch_id, "planned_status": "planned", "data_path": str(branch_source),
                        "source_sha256": branch_source_sha256, "config_path": str(config_path),
                        "config_file_sha256": study_service._hash_file(config_path),
                        "config_sha256": clustering_config_sha256(branch_config),
                        "scientific_config_sha256": clustering_scientific_config_sha256(branch_config),
                        "manifest_path": str(manifest_path), "manifest_file_sha256": study_service._hash_file(manifest_path),
                        "manifest_sha256": stable_fingerprint(branch_manifest.model_dump(mode="json")),
                        "run_dir": str(root / "runs" / branch_id),
                        "selected_ids_sha256": spec["selected_ids_sha256"], "raw_counts_sha256": spec["raw_counts_sha256"],
                        "feature_axis_sha256": spec["feature_axis_sha256"],
                        "raw_role_columns": role_columns, "raw_roles_sha256": stable_fingerprint(role_values),
                        "predicted_retained_ids_sha256": spec["predicted_retained_ids_sha256"]})
    del canonical
    if study_service._hash_file(source) != source_sha256 or execution_environment() != environment:
        raise ValueError("Original source or execution environment changed while preparing frozen inputs")
    bindings = {"schema_version": 1, "policy_version": plan["policy_version"], "output_dir": str(root),
                "plan_path": str(root / "plan.json"), "plan_file_sha256": study_service._hash_file(root / "plan.json"),
                "plan_sha256": plan["plan_sha256"], "source_path": str(source), "source_sha256": source_sha256,
                "execution_environment": environment, "branches": records,
                "preparation_scope": "Raw source selection and derived-file roundtrip verification only; no native computation or reference evaluation",
                "derived_lineage_scope": "File SHA, complete ordered cell mask, exact count/dtype and feature-axis identity; no cross-registry parent URI is invented"}
    _write(root / "bindings.json", bindings)
    return {"output_dir": str(root), "bindings_path": str(root / "bindings.json"),
            "bindings_sha256": study_service._hash_file(root / "bindings.json"),
            "plan_path": str(root / "plan.json"), "plan_file_sha256": bindings["plan_file_sha256"]}


def _load_preparation(prepared, expected_bindings_sha256=None):
    from eacbp.application import study_service

    if isinstance(prepared, Mapping):
        path = Path(prepared["bindings_path"]).expanduser().resolve()
        expected = prepared["bindings_sha256"]
        if expected_bindings_sha256 is not None and expected != expected_bindings_sha256:
            raise ValueError("Clustering preparation SHA differs from the supplied protocol pin")
    else:
        path = Path(prepared).expanduser().resolve()
        path = path / "bindings.json" if path.is_dir() else path
        expected = expected_bindings_sha256
        if expected is None:
            raise ValueError("A saved preparation requires its pre-execution bindings SHA")
    if study_service._hash_file(path) != expected:
        raise ValueError("Frozen clustering bindings changed before execution")
    bindings = _read(path)
    root = Path(bindings["output_dir"]).resolve()
    if path != root / "bindings.json":
        raise ValueError("Clustering bindings must belong to their frozen output directory")
    plan_path = _contained(bindings["plan_path"], root)
    if study_service._hash_file(plan_path) != bindings["plan_file_sha256"]:
        raise ValueError("Frozen clustering plan file changed")
    plan = _read(plan_path)
    validate_clustering_sensitivity_plan(plan)
    if plan["plan_sha256"] != bindings["plan_sha256"] or plan["source_sha256"] != bindings["source_sha256"]:
        raise ValueError("Clustering preparation and plan disagree")
    expected_ids = [item["branch_id"] for item in [plan["primary"], *plan["branches"]]]
    if [item["branch_id"] for item in bindings["branches"]] != expected_ids:
        raise ValueError("Frozen preparation omitted, added or reordered branches")
    return root, bindings, plan, expected


def _validate_frozen_inputs(root, bindings, plan):
    from eacbp.application import study_service

    if execution_environment() != bindings["execution_environment"]:
        raise ValueError("Execution environment differs from pre-native preparation")
    if study_service._hash_file(Path(bindings["source_path"])) != bindings["source_sha256"]:
        raise ValueError("Original raw source changed after clustering preparation")
    verified_sources = {}
    for record in bindings["branches"]:
        spec = _spec(plan, record["branch_id"])
        if record["planned_status"] != spec["planned_status"]:
            raise ValueError("Frozen branch execution eligibility changed")
        if record["planned_status"] != "planned":
            continue
        for name in ("config", "manifest"):
            path = _contained(record[name + "_path"], root)
            if study_service._hash_file(path) != record[name + "_file_sha256"]:
                raise ValueError(f"Frozen branch {name} file changed")
        config = _read(record["config_path"])
        validate_clustering_branch_config(plan, record["branch_id"], config)
        if clustering_config_sha256(config) != record["config_sha256"]:
            raise ValueError("Frozen branch configuration binding changed")
        saved_manifest = StudyManifest.model_validate(_read(record["manifest_path"]))
        if saved_manifest.data.raw_artifact_uri is not None or stable_fingerprint(saved_manifest.model_dump(mode="json")) != record["manifest_sha256"]:
            raise ValueError("Frozen branch manifest identity changed")
        _contained(record["run_dir"], root)
        source = Path(record["data_path"]).resolve()
        if spec["source_kind"] == "original":
            if source != Path(bindings["source_path"]).resolve() or record["source_sha256"] != bindings["source_sha256"]:
                raise ValueError("Full-source branch is not bound to the original source")
        else:
            _contained(source, root / "inputs")
        if source not in verified_sources:
            verified_sources[source] = study_service._hash_file(source)
        if verified_sources[source] != record["source_sha256"]:
            raise ValueError("Frozen native input source file changed")
        if config["autonomous_analysis"]["source_sha256"] != record["source_sha256"]:
            raise ValueError("Frozen native configuration is bound to another input source")


def _native_file_hashes(path):
    from eacbp.application import study_service

    return {str(item.resolve()): study_service._hash_file(item)
            for pattern in ("*.json", "report.md", "snapshots/*.json", "summaries/*.json", "*.jsonl", "artifacts/.artifact_registry.json")
            for item in path.glob(pattern) if item.is_file()}


def _verify_upstream_metadata(task, contract, metadata, payload):
    parameters = contract.get("parameters", {})
    if task.capability == "normalization":
        expected = {"target_sum": float(parameters.get("target_sum", 10000.)),
                    "n_top_genes": int(parameters.get("n_top_genes", 300)),
                    "input_layer": parameters.get("input_layer", "counts")}
        actual = {name: metadata.parameters.get(name) for name in expected}
        normalization = model_to_jsonable(payload.uns.get("normalization", {}))
        if (actual != expected or metadata.parameters.get("counts_layer") != "counts"
                or normalization.get("method") != task.method_used
                or normalization.get("target_sum") != expected["target_sum"]
                or normalization.get("input_layer") != expected["input_layer"]
                or normalization.get("log1p") is not True
                or normalization.get("counts_layer") != "counts"):
            raise ValueError("Native normalization metadata/receipt differ from the durable contract")
    elif task.capability == "integration" and task.method_used == "no_correction_v1":
        expected = {"n_components": int(parameters.get("n_components", 20)),
                    "random_seed": int(parameters.get("random_seed", 42))}
        representation = payload.uns.get("integration", {}).get("representation_receipt", {})
        if ({name: metadata.parameters.get(name) for name in expected} != expected
                or representation.get("requested_n_components") != expected["n_components"]
                or representation.get("random_seed") != expected["random_seed"]):
            raise ValueError("Native integration metadata/representation receipt differ from the durable contract")


def collect_clustering_run(run_dir, *, plan, binding) -> dict[str, Any]:
    """Read every output through durable admission before supplying comparison labels."""
    from eacbp.application import study_service
    from eacbp.auditor.clustering import validate_clustering_execution_receipt
    from eacbp.capabilities.sc_data import SCData
    from eacbp.evidence.snapshot import load_study_snapshot

    path = Path(run_dir).resolve()
    spec = _spec(plan, binding["branch_id"])
    registry = study_service._make_registry(path / "artifacts")
    snapshot = load_study_snapshot(path / "snapshot.json", artifact_registry=registry)
    if clustering_config_sha256(snapshot.config) != binding["config_sha256"]:
        raise ValueError("Native snapshot configuration differs from pre-execution binding")
    validate_clustering_branch_config(plan, spec["branch_id"], snapshot.config)
    expected_manifest = _read(binding["manifest_path"])
    actual_manifest = snapshot.manifest.model_dump(mode="json")
    actual_manifest["data"]["raw_artifact_uri"] = None
    if stable_fingerprint(actual_manifest) != stable_fingerprint(expected_manifest):
        raise ValueError("Native snapshot manifest differs from its frozen study")
    provenance = {"run_dir": str(path), "report": str(path / "report.md"),
                  "snapshot": str(path / "snapshot.json"), "files": _native_file_hashes(path),
                  "artifacts": [], "audit_signatures": []}
    base = {"admitted": False, "run_status": snapshot.summary.get("status"),
            "analysis_request_status": snapshot.summary.get("analysis_request_status"), "provenance": provenance}
    if snapshot.summary.get("status") != "success" or snapshot.summary.get("analysis_request_status") != "complete":
        return {**base, "reason": "Native execution or source-bound scientific request is incomplete"}
    required = {"dataset_audit", "qc", "qc_review", "normalization", "integration", "clustering", "cluster_markers"}
    counts = {name: sum(task.capability == name for task in snapshot.task_results) for name in required}
    if set(task.capability for task in snapshot.task_results) != required or any(value != 1 for value in counts.values()):
        raise ValueError("Native clustering replay lacks exactly one complete task per frozen capability")
    records = registry.list_audit_records()
    chain_data = {}
    clustering_execution_receipt = None
    for task in snapshot.task_results:
        if task.status.value != "success":
            raise ValueError("A native clustering task did not succeed")
        result_identity = registry._result_identity(task.model_dump(mode="python"))
        matches = [record for record in records if record.task_id == task.task_id and record.overall_passed
                   and not record.stop_rule_triggered and record.result_fingerprint == result_identity
                   and record.output_artifacts == task.output_artifacts]
        if len(matches) != 1:
            raise ValueError("Native task lacks unique matching durable audit admission")
        record = matches[0]
        for name, value in snapshot.config.get("capability_parameters", {}).get(task.capability, {}).items():
            if stable_fingerprint(record.contract.get("parameters", {}).get(name)) != stable_fingerprint(value):
                raise ValueError("Actual native task parameters differ from the frozen configuration")
        outputs = registry.get_audited_many(record.signature, record.contract)
        provenance["audit_signatures"].append(record.signature)
        for metadata, payload in outputs:
            provenance["artifacts"].append({"uri": metadata.uri, "sha256": metadata.sha256_hash,
                                            "path": str(registry.storage.resolve_path(metadata.storage_path)),
                                            "task_id": task.task_id, "audit_signature": record.signature})
            if task.capability in {"qc", "normalization", "integration", "clustering"}:
                scientific_data = SCData.from_payload(payload)
                _verify_upstream_metadata(task, record.contract, metadata, scientific_data)
                chain_data[task.capability] = scientific_data
        if task.capability == "clustering":
            clustering_execution_receipt = validate_clustering_execution_receipt(
                TaskContract.model_validate(record.contract), task, registry, recompute_partition=False)
    raw_uri = snapshot.manifest.data.raw_artifact_uri
    raw_metadata, raw = registry.get_scdata(raw_uri)
    if (raw_metadata.parameters.get("source_sha256") != binding["source_sha256"]
            or Path(raw_metadata.parameters.get("source_file", "")).resolve() != Path(binding["data_path"]).resolve()):
        raise ValueError("Native ingest source differs from pre-execution raw-file identity")
    _verify_subset(raw, spec)
    role_columns = binding["raw_role_columns"]
    if not set(role_columns) <= set(raw.obs) or stable_fingerprint(_role_values(raw.obs, role_columns)) != binding["raw_roles_sha256"]:
        raise ValueError("Native raw biological/technical metadata differ from frozen prepared source roles")
    qc = chain_data["qc"]
    retained_ids = qc.obs.index.astype(str).tolist()
    if retained_ids != spec["predicted_retained_ids"]:
        raise ValueError("Native QC retention differs from frozen full-source thresholds")
    actual_qc = {key: qc.uns["qc"].get(key) for key in ("min_genes", "max_mito_pct")}
    if actual_qc != spec["qc_parameters"]:
        raise ValueError("Native QC used changed thresholds")
    positions = raw.obs.index.astype(str).get_indexer(retained_ids)
    expected_counts = raw_count_fingerprint(_counts(raw)[positions])
    expected_roles = _role_values(raw.obs.iloc[positions], role_columns)
    for capability, payload in chain_data.items():
        if payload.obs.index.astype(str).tolist() != retained_ids or raw_count_fingerprint(_counts(payload)) != expected_counts:
            raise ValueError(f"Native {capability} changed admitted raw counts or retained cell order")
        if feature_axis_fingerprint(payload.var, "var.gene_name") != spec["feature_axis_sha256"]:
            raise ValueError(f"Native {capability} changed the original complete feature universe")
        if not set(role_columns) <= set(payload.obs) or stable_fingerprint(_role_values(payload.obs, role_columns)) != stable_fingerprint(expected_roles):
            raise ValueError(f"Native {capability} changed admitted source role values")
    clusters = chain_data["clustering"].obs["cluster"].copy()
    actual_clustering = {"resolution": clustering_execution_receipt["leiden"]["resolution"],
                         "random_seed": clustering_execution_receipt["leiden"]["random_state"],
                         "n_neighbors": clustering_execution_receipt["neighbors"]["requested_n_neighbors"]}
    if actual_clustering != spec["clustering_parameters"]:
        raise ValueError("Actual native clustering receipt differs from the frozen perturbation")
    return {**base, "admitted": True, "source_sha256": binding["source_sha256"],
            "expected_source_sha256": binding["source_sha256"], "config_sha256": clustering_config_sha256(snapshot.config),
            "expected_config_sha256": binding["config_sha256"],
            "scientific_config_sha256": clustering_scientific_config_sha256(snapshot.config),
            "raw_selected_ids": raw.obs.index.astype(str).tolist(), "selected_ids_sha256": spec["selected_ids_sha256"],
            "raw_counts_sha256": raw_count_fingerprint(_counts(raw)),
            "feature_axis_sha256": feature_axis_fingerprint(raw.var, "var.gene_name"),
            "qc_parameters": actual_qc, "clustering_parameters": actual_clustering,
            "retained_ids": retained_ids, "clusters": clusters,
            "clustering_execution_receipt": clustering_execution_receipt}


def _markdown(result):
    lines = ["# Frozen clustering sensitivity", "", f"Execution/comparison status: `{result['status']}`.", "",
             result.get("interpretation", "Native failures and missing branches remain explicit."), "",
             "| Branch | Native status | Comparison | Shared cells | ARI | NMI |", "| :--- | :--- | :--- | ---: | ---: | ---: |"]
    for item in result.get("branches", []):
        lines.append(f"| {item['branch_id']} | {item.get('run_status', 'not_run')} | {item['status']} | {item.get('n_shared_cells', 'unavailable')} | {item.get('shared_cell_ari', 'unavailable')} | {item.get('shared_cell_nmi', 'unavailable')} |")
    lines.extend(["", "All original primary clusters, matching denominators, missing comparisons, native audit receipts and artifact SHA values are retained in comparison.json.",
                  "Completion describes bounded native execution and paired comparisons. Biological validation and any separately frozen scientific acceptance are distinct.", ""])
    return "\n".join(lines)


def execute_clustering_sensitivity(prepared, *, expected_bindings_sha256=None, resume=False) -> dict[str, Any]:
    """Execute pinned studies; explicit resume re-admits completed native outputs.

    Resume never treats journal status or a saved comparison as scientific
    admission. Incomplete existing studies remain preserved and unavailable;
    task-level recovery stays with study_service.resume_study.
    """
    from eacbp.application.execution_journal import ExecutionJournal
    from eacbp.application.process_resources import ProcessTreeSampler

    if type(resume) is not bool:
        raise ValueError("Clustering resume must be an explicit boolean")
    root, bindings, plan, bindings_sha256 = _load_preparation(prepared, expected_bindings_sha256)
    if not resume and any((root / name).exists() for name in ("comparison.json", "runs", "execution_journal")):
        raise ValueError("Frozen clustering preparation was already executed; replays require a fresh directory")
    identity = {key: bindings[key] for key in ("plan_file_sha256", "plan_sha256", "source_sha256", "execution_environment")}
    identity["bindings_sha256"] = bindings_sha256
    with ExecutionJournal(root, identity=identity, resume=resume) as journal:
        invocation = journal.start_attempt("invocation", details={"resume": resume})
        resource_dir = journal.path / "resources" / invocation
        try:
            with ProcessTreeSampler(resource_dir) as sampler:
                output = _execute_clustering_stages(root, bindings, plan, bindings_sha256,
                                                     journal=journal, invocation=invocation, resume=resume)
                journal.finish_attempt(invocation, status="completed" if output["status"] == "completed" else "failed",
                                       details={"comparison": output["comparison"], "status": output["status"]})
        except BaseException as exc:
            if next(attempt for attempt in journal.attempts if attempt["attempt_id"] == invocation)["status"] == "started":
                journal.finish_attempt(invocation, status="failed", details={"error_type": type(exc).__name__, "reason": str(exc)})
            raise
        output.update(execution_journal=str(journal.path), execution_journal_identity_sha256=journal.identity_sha256,
                      invocation_id=invocation, process_resources=str(sampler.summary_path) if sampler.summary_path.is_file() else None,
                      process_resource_samples=str(sampler.sample_path),
                      process_resource_summary_persisted=sampler.summary_path.is_file())
        return output


def _execute_clustering_stages(root, bindings, plan, bindings_sha256, *, journal, invocation, resume):
    from eacbp.application import study_service
    from eacbp.evidence.snapshot import load_study_snapshot

    receipts, failures, reused, started, prior_pins = {}, [], [], [], {}
    preflight = journal.start_attempt("preflight")
    try:
        _validate_frozen_inputs(root, bindings, plan)
        for binding in bindings["branches"]:
            branch_id = binding["branch_id"]
            if binding["planned_status"] != "planned":
                continue
            previous = [attempt for attempt in journal.attempts if attempt["branch_id"] == branch_id
                        and attempt["stage"] in {"native_execution", "native_admission"}
                        and attempt["status"] == "completed"
                        and attempt["finish"]["details"].get("admitted") is True]
            if previous and not Path(binding["run_dir"]).exists():
                raise ValueError("Previously admitted native directory is missing; resume cannot recompute it")
            if not Path(binding["run_dir"]).exists():
                continue
            prior_pins[branch_id] = previous[-1]["finish"]["details"]["provenance"]["files"] if previous else None
            if prior_pins[branch_id] is not None and _native_file_hashes(Path(binding["run_dir"])) != prior_pins[branch_id]:
                raise ValueError("Previously admitted native evidence changed before resume")
        journal.finish_attempt(preflight, status="completed")
        for binding in bindings["branches"]:
            branch_id = binding["branch_id"]
            if binding["planned_status"] != "planned":
                continue
            if branch_id != "primary" and not receipts.get("primary", {}).get("admitted"):
                failures.append({"branch_id": branch_id, "status": "not_run", "reason": "Primary native study was not fully admitted; no alternate branch rescues it"})
                continue
            _validate_frozen_inputs(root, bindings, plan)
            config = _read(binding["config_path"])
            manifest = StudyManifest.model_validate(_read(binding["manifest_path"]))
            branch_dir = Path(binding["run_dir"])
            existing = branch_dir.exists()
            if prior_pins.get(branch_id) is not None and _native_file_hashes(branch_dir) != prior_pins[branch_id]:
                raise ValueError("Previously admitted native evidence changed before collection")
            native_attempt = journal.start_attempt("native_admission" if existing else "native_execution", branch_id)
            try:
                if existing:
                    if not resume:
                        raise ValueError("Existing native output requires explicit resume")
                else:
                    started.append(branch_id)
                    study_service.run_study(manifest=manifest, config=config, data=binding["data_path"], run_dir=branch_dir)
                receipts[branch_id] = collect_clustering_run(branch_dir, plan=plan, binding=binding)
                if existing and receipts[branch_id]["admitted"] is True:
                    reused.append(branch_id)
                if receipts[branch_id]["admitted"] is not True:
                    failures.append({"branch_id": branch_id, "status": "incomplete", "reason": receipts[branch_id]["reason"], "run_dir": str(branch_dir)})
                journal.finish_attempt(native_attempt, status="completed" if receipts[branch_id]["admitted"] is True else "failed",
                                       details={"admitted": receipts[branch_id]["admitted"], "run_dir": str(branch_dir),
                                                "reused": branch_id in reused, "provenance": receipts[branch_id]["provenance"]})
            except Exception as exc:
                receipts[branch_id] = {"admitted": False, "run_status": "failed",
                                       "reason": f"{type(exc).__name__}: {exc}",
                                       "provenance": {"run_dir": str(branch_dir), "files": _native_file_hashes(branch_dir)}}
                failures.append({"branch_id": branch_id, "status": "failed", "reason": receipts[branch_id]["reason"], "run_dir": str(branch_dir)})
                journal.finish_attempt(native_attempt, status="failed", details={"reason": receipts[branch_id]["reason"],
                                       "run_dir": str(branch_dir), "reused": branch_id in reused})
        comparison_attempt = journal.start_attempt("comparison")
        _validate_frozen_inputs(root, bindings, plan)
        if study_service._hash_file(root / "bindings.json") != bindings_sha256 or study_service._hash_file(root / "plan.json") != bindings["plan_file_sha256"]:
            raise ValueError("Frozen preparation changed during native execution")
        for receipt in receipts.values():
            provenance = receipt["provenance"]
            for path, digest in provenance.get("files", {}).items():
                if study_service._hash_file(Path(path)) != digest:
                    raise ValueError("Native primary or branch evidence changed during sensitivity execution")
            if receipt["admitted"] is True:
                native = Path(provenance["run_dir"])
                load_study_snapshot(native / "snapshot.json", artifact_registry=study_service._make_registry(native / "artifacts"))
        result = compare_clustering_sensitivity_runs(plan, baseline=receipts.get("primary", {"admitted": False}),
                                                     branches={key: value for key, value in receipts.items() if key != "primary"})
        result.update(failures=failures, bindings_sha256=bindings_sha256, plan_file_sha256=bindings["plan_file_sha256"],
                      primary_native_summary={key: value for key, value in receipts.get("primary", {}).items()
                                              if key in {"run_status", "analysis_request_status", "admitted", "reason"}},
                      native_execution_receipts={key: value.get("clustering_execution_receipt") for key, value in receipts.items()},
                      run_provenance={key: value["provenance"] for key, value in receipts.items()},
                      native_evidence_preserved=True, original_source_preserved=True,
                      execution_environment_preserved=True, references_opened=False,
                      resume_prior_native_file_pins={key: value is not None for key, value in prior_pins.items()},
                      resumed_without_previous_native_file_pin=[key for key in reused if prior_pins.get(key) is None],
                      resume_integrity_scope="Re-admission is mandatory; journal-pinned prior files are exact. A crash before first collection has no prior sensitivity file pin and relies on durable native task evidence plus the frozen inputs.")
    except Exception as exc:
        for attempt in journal.attempts:
            if attempt["status"] == "started" and attempt["stage"] in {"preflight", "comparison"}:
                journal.finish_attempt(attempt["attempt_id"], status="failed", details={"reason": f"{type(exc).__name__}: {exc}"})
        result = {"status": "failed", "plan_id": plan["plan_id"], "plan_sha256": plan["plan_sha256"],
                  "bindings_sha256": bindings_sha256, "reason": f"{type(exc).__name__}: {exc}",
                  "failures": failures, "run_provenance": {key: value["provenance"] for key, value in receipts.items()},
                  "references_opened": False, "winner_selected": False, "primary_replaced": False,
                  "biological_stability_claim_supported": False}
    output_dir = root / "comparisons" / invocation
    comparison_path = output_dir / "comparison.json"
    try:
        output_dir.mkdir(parents=True, exist_ok=False)
        _write(comparison_path, result)
        if not (root / "comparison.json").exists():
            _write(root / "comparison.json", result)
        for attempt in journal.attempts:
            if attempt["status"] == "started" and attempt["stage"] == "comparison":
                journal.finish_attempt(attempt["attempt_id"], status="completed", details={"status": result["status"],
                                       "path": str(comparison_path), "sha256": study_service._hash_file(comparison_path)})
    except BaseException as exc:
        for attempt in journal.attempts:
            if attempt["status"] == "started" and attempt["stage"] == "comparison":
                journal.finish_attempt(attempt["attempt_id"], status="failed",
                                       details={"error_type": type(exc).__name__, "reason": str(exc),
                                                "path": str(comparison_path)})
        raise
    postprocessing = journal.start_attempt("postprocessing")
    report_path = output_dir / "report.md"
    try:
        report_path.write_text(_markdown(result), encoding="utf-8")
        if not (root / "report.md").exists():
            (root / "report.md").write_text(report_path.read_text(encoding="utf-8"), encoding="utf-8")
        _write(root / "latest_comparison.json", {"invocation_id": invocation, "comparison": str(comparison_path),
               "comparison_sha256": study_service._hash_file(comparison_path), "report": str(report_path),
               "report_sha256": study_service._hash_file(report_path), "status": result["status"]})
        journal.finish_attempt(postprocessing, status="completed", details={"report": str(report_path),
                               "sha256": study_service._hash_file(report_path)})
    except BaseException as exc:
        journal.finish_attempt(postprocessing, status="postprocessing_failed",
                               details={"error_type": type(exc).__name__, "reason": str(exc), "comparison": str(comparison_path)})
        raise
    return {"status": result["status"], "output_dir": str(root), "plan_id": plan["plan_id"],
            "bindings_sha256": bindings_sha256, "comparison": str(comparison_path), "report": str(report_path),
            "native_reused_branch_ids": reused, "native_started_branch_ids": started,
            "primary_native_summary": result.get("primary_native_summary"),
            "completion_scope": "Detached prespecified native execution and paired descriptions only",
            "biological_stability_claim_supported": False}


def run_clustering_sensitivity(*, manifest, config, data, output_dir):
    """Prepare and execute one bounded standalone experiment through native run_study."""
    prepared = prepare_clustering_sensitivity(manifest=manifest, config=config, data=data, output_dir=output_dir)
    return execute_clustering_sensitivity(prepared)
