"""Frozen input tampering, failure isolation, and genuine native admission checks."""

from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.application import clustering_sensitivity_execution as execution
from eacbp.application import study_service
from eacbp.application.analysis_requirements import freeze_question_requirements
from eacbp.application.autonomous import observe_input
from eacbp.application.clustering_sensitivity import (
    clustering_config_sha256, clustering_scientific_config_sha256,
)
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.study import BiologicalDesign, ExperimentalDesign, StudyManifest


def _source(tmp_path, *, donors=0):
    rng = np.random.default_rng(35)
    values = rng.poisson(1, (120, 40)).astype(np.int32)
    for group in range(4):
        values[group * 30:(group + 1) * 30, group * 8:(group + 1) * 8] += rng.poisson(15, (30, 8))
    values[:3] = 0
    counts = sparse.csr_matrix(values)
    obs = pd.DataFrame(index=[f"cell_{number}" for number in range(120)])
    if donors:
        obs["donor_id"] = [f"GSM{2230757 + number % donors}" for number in range(120)]
    var = pd.DataFrame({"gene_name": [f"GENE{number}" for number in range(40)]},
                       index=[f"feature{number}" for number in range(40)])
    payload = SCData(counts, obs, var, layers={"counts": counts})
    path = tmp_path / "raw-source.h5ad"
    payload.to_anndata().write_h5ad(path)
    observation = observe_input(payload, question="Explore cell populations and marker genes")
    source_sha256 = study_service._hash_file(path)
    settings = {"intent": "exploration", "question": "Explore cell populations and marker genes",
                "source_sha256": source_sha256, "observations": observation,
                "policy_version": "observed_single_cell_v3",
                "input_preparation": {key: observation[key] for key in ("counts_source", "gene_name_source", "obs_column_map")},
                "requirement_contract": freeze_question_requirements(question="Explore cell populations and marker genes",
                                                                       observations=observation, source_sha256=source_sha256)}
    config = {"mode": "real", "method_profile": "standard", "threads": 1,
              "autonomous_analysis": settings, "analysis_extensions": {"cluster_markers": True, "qc_review": True},
              "method_overrides": {"integration": "no_correction_v1"},
              "capability_parameters": {"qc": {key: observation["qc"][key] for key in ("min_genes", "max_mito_pct")},
                                        "normalization": {"n_top_genes": 35, "target_sum": 10000.},
                                        "integration": {"n_components": 10, "random_seed": 42},
                                        "clustering": {"resolution": 1., "random_seed": 42, "n_neighbors": 15}}}
    manifest = StudyManifest(study_id="source_replay", biological_design=BiologicalDesign(species="human", tissue="test"),
                             experimental_design=ExperimentalDesign(biological_unit="donor" if donors else "unknown",
                                                                    total_samples=donors))
    return path, config, manifest


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    monkeypatch.setattr(execution, "execution_environment", lambda: {"source": "frozen-test-environment", "source_sha256": "e" * 64})
    source, config, manifest = _source(tmp_path)
    record = execution.prepare_clustering_sensitivity(manifest=manifest, config=config, data=source,
                                                     output_dir=tmp_path / "experiment")
    return record, source, config, manifest


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _fake_receipt(plan, binding, path):
    spec = next(item for item in [plan["primary"], *plan["branches"]] if item["branch_id"] == binding["branch_id"])
    config = _read(binding["config_path"])
    ids = spec["predicted_retained_ids"]
    labels = pd.Series([str(int(cell.split("_")[1]) % 3) for cell in ids], index=ids)
    return {"admitted": True, "run_status": "success", "analysis_request_status": "complete",
            "source_sha256": binding["source_sha256"], "expected_source_sha256": binding["source_sha256"],
            "config_sha256": clustering_config_sha256(config), "expected_config_sha256": binding["config_sha256"],
            "scientific_config_sha256": clustering_scientific_config_sha256(config),
            "raw_selected_ids": spec["selected_ids"], "selected_ids_sha256": spec["selected_ids_sha256"],
            "raw_counts_sha256": spec["raw_counts_sha256"], "feature_axis_sha256": spec["feature_axis_sha256"],
            "qc_parameters": spec["qc_parameters"], "clustering_parameters": spec["clustering_parameters"],
            "retained_ids": ids, "clusters": labels,
            "provenance": {"run_dir": str(path), "files": execution._native_file_hashes(path)}}


def _stub_native(monkeypatch, *, hook=None, failed=None, incomplete=None):
    calls = []

    def run(**kwargs):
        branch_id = kwargs["config"]["autonomous_analysis"]["clustering_sensitivity_branch"]["branch_id"]
        calls.append(branch_id)
        path = Path(kwargs["run_dir"])
        path.mkdir(parents=True)
        (path / "summary.json").write_text(json.dumps({"branch_id": branch_id, "status": "success"}), encoding="utf-8")
        if hook:
            hook(branch_id, kwargs)
        if branch_id == failed:
            raise RuntimeError("Preserved native failure")
        return {"status": "success"}

    def collect(path, *, plan, binding):
        receipt = _fake_receipt(plan, binding, path)
        if binding["branch_id"] == incomplete:
            receipt.update(admitted=False, reason="Native requested capability is incomplete")
        return receipt

    monkeypatch.setattr(study_service, "run_study", run)
    monkeypatch.setattr(execution, "collect_clustering_run", collect)
    monkeypatch.setattr("eacbp.evidence.snapshot.load_study_snapshot", lambda *args, **kwargs: None)
    monkeypatch.setattr(study_service, "_make_registry", lambda path: None)
    return calls


def test_preparation_writes_and_roundtrips_every_source_before_native(tmp_path, monkeypatch):
    monkeypatch.setattr(execution, "execution_environment", lambda: {"source": "frozen", "source_sha256": "e" * 64})
    source, config, manifest = _source(tmp_path, donors=4)
    source_before = study_service._hash_file(source)
    monkeypatch.setattr(study_service, "run_study", lambda **kwargs: pytest.fail("Preparation must never execute a native study"))
    prepared = execution.prepare_clustering_sensitivity(manifest=manifest, config=config, data=source,
                                                        output_dir=tmp_path / "experiment")
    bindings = _read(prepared["bindings_path"])
    assert len(bindings["branches"]) == 9
    assert source_before == study_service._hash_file(source)
    assert not (Path(prepared["output_dir"]) / "runs").exists()
    assert len(list((Path(prepared["output_dir"]) / "inputs").glob("*.h5ad"))) == 5
    for binding in bindings["branches"]:
        saved_manifest = _read(binding["manifest_path"])
        assert saved_manifest["data"]["raw_artifact_uri"] is None
        if binding["branch_id"].startswith("leave_out"):
            assert saved_manifest["experimental_design"]["total_samples"] == 3
        assert {"donor_id", "donor"} <= set(binding["raw_role_columns"])
        assert study_service._hash_file(Path(binding["data_path"])) == binding["source_sha256"]
    assert len({_read(item["manifest_path"])["study_id"] for item in bindings["branches"]}) == 9


@pytest.mark.parametrize("target", ["plan", "bindings", "config", "manifest", "derived_source", "original_source"])
def test_mutation_before_first_native_is_rejected_without_any_execution(prepared, monkeypatch, target):
    record, source, _, _ = prepared
    bindings = _read(record["bindings_path"])
    if target == "plan":
        path = Path(record["plan_path"])
    elif target == "bindings":
        path = Path(record["bindings_path"])
    elif target in {"config", "manifest"}:
        path = Path(bindings["branches"][1][target + "_path"])
    elif target == "derived_source":
        path = Path(bindings["branches"][-1]["data_path"])
    else:
        path = source
    path.write_bytes(path.read_bytes() + b"changed")
    calls = _stub_native(monkeypatch)
    if target in {"plan", "bindings"}:
        with pytest.raises(ValueError):
            execution.execute_clustering_sensitivity(record)
    else:
        result = execution.execute_clustering_sensitivity(record)
        assert result["status"] == "failed"
    assert calls == []


def test_environment_change_before_first_native_is_explicit_failure(prepared, monkeypatch):
    record, _, _, _ = prepared
    calls = _stub_native(monkeypatch)
    monkeypatch.setattr(execution, "execution_environment", lambda: {"source": "changed", "source_sha256": "f" * 64})
    result = execution.execute_clustering_sensitivity(record)
    assert result["status"] == "failed" and calls == []
    assert "environment" in _read(result["comparison"])["reason"]


def test_direct_native_replays_are_bounded_and_keep_primary_and_qc_unchanged(prepared, monkeypatch):
    record, source, config, _ = prepared
    source_before = study_service._hash_file(source)
    config_before = deepcopy(config)
    calls = _stub_native(monkeypatch)
    result = execution.execute_clustering_sensitivity(record)
    assert calls == ["primary", "resolution_0_5", "resolution_2_0", "clustering_seed_17", "uniform_sample_80pct"]
    assert result["status"] == "completed"
    comparison = _read(result["comparison"])
    assert comparison["primary_replaced"] is False and comparison["winner_selected"] is False
    assert comparison["references_opened"] is False
    assert comparison["execution_environment_preserved"] is True
    assert source_before == study_service._hash_file(source) and config == config_before
    with pytest.raises(ValueError, match="already executed"):
        execution.execute_clustering_sensitivity(record)


@pytest.mark.parametrize("mode", ["failure", "incomplete"])
def test_failed_or_incomplete_primary_is_preserved_and_never_rescued(prepared, monkeypatch, mode):
    record, _, _, _ = prepared
    calls = _stub_native(monkeypatch, failed="primary" if mode == "failure" else None,
                         incomplete="primary" if mode == "incomplete" else None)
    result = execution.execute_clustering_sensitivity(record)
    assert calls == ["primary"] and result["status"] == "partial"
    comparison = _read(result["comparison"])
    assert len(comparison["failures"]) == 5
    assert comparison["primary_native_summary"]["admitted"] is False
    assert all(item["status"] == "not_available" for item in comparison["branches"])
    assert (Path(record["output_dir"]) / "runs" / "primary" / "summary.json").is_file()


def test_one_failed_branch_does_not_remove_it_or_select_another_winner(prepared, monkeypatch):
    record, _, _, _ = prepared
    calls = _stub_native(monkeypatch, failed="resolution_0_5")
    result = execution.execute_clustering_sensitivity(record)
    assert len(calls) == 5 and result["status"] == "partial"
    comparison = _read(result["comparison"])
    failed = next(item for item in comparison["branches"] if item["branch_id"] == "resolution_0_5")
    assert failed["status"] == "not_available" and failed["run_status"] == "failed"
    assert comparison["winner_selected"] is False
    assert comparison["failures"][0]["branch_id"] == "resolution_0_5"


def test_later_branch_cannot_rewrite_primary_evidence(prepared, monkeypatch):
    record, _, _, _ = prepared
    primary = Path(record["output_dir"]) / "runs" / "primary" / "summary.json"

    def hook(branch_id, kwargs):
        if branch_id == "resolution_0_5":
            primary.write_text("changed primary", encoding="utf-8")

    _stub_native(monkeypatch, hook=hook)
    result = execution.execute_clustering_sensitivity(record)
    assert result["status"] == "failed"
    assert "Native primary or branch evidence changed" in _read(result["comparison"])["reason"]


def test_saved_preparation_requires_the_original_pinned_digest(prepared, monkeypatch):
    record, _, _, _ = prepared
    calls = _stub_native(monkeypatch)
    with pytest.raises(ValueError, match="pre-execution bindings SHA"):
        execution.execute_clustering_sensitivity(record["output_dir"])
    result = execution.execute_clustering_sensitivity(record["output_dir"], expected_bindings_sha256=record["bindings_sha256"])
    assert result["status"] == "completed" and len(calls) == 5


def test_genuine_native_primary_is_collectible_only_with_full_durable_admission(prepared, monkeypatch):
    record, source, _, _ = prepared
    bindings = _read(record["bindings_path"])
    plan = _read(record["plan_path"])
    primary = bindings["branches"][0]
    manifest = StudyManifest.model_validate(_read(primary["manifest_path"]))
    source_before = study_service._hash_file(source)
    summary = study_service.run_study(manifest=manifest, config=_read(primary["config_path"]),
                                      data=primary["data_path"], run_dir=primary["run_dir"])
    assert summary["status"] == "success" and summary["analysis_request_status"] == "complete"
    receipt = execution.collect_clustering_run(primary["run_dir"], plan=plan, binding=primary)
    assert receipt["admitted"] is True
    assert len(receipt["provenance"]["audit_signatures"]) == 7
    assert len(receipt["provenance"]["artifacts"]) == 10
    assert receipt["retained_ids"] == plan["primary"]["predicted_retained_ids"]
    assert receipt["clustering_execution_receipt"]["leiden"]["random_state"] == 42
    assert study_service._hash_file(source) == source_before
    changed = deepcopy(primary)
    changed["config_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="configuration"):
        execution.collect_clustering_run(primary["run_dir"], plan=plan, binding=changed)
    # Every sibling is genuinely admitted and its payload byte hash remains
    # valid. Returned metadata cannot override the durable scientific contract.
    registry = study_service._make_registry(Path(primary["run_dir"]) / "artifacts")
    for capability, key, value in (("normalization", "target_sum", 20000.),
                                   ("integration", "random_seed", 17)):
        admitted_reader = registry.get_audited_many

        def altered_metadata(signature, contract):
            outputs = admitted_reader(signature, contract)
            if contract["capability"] == capability:
                for metadata, _ in outputs:
                    metadata.parameters[key] = value
            return outputs

        with monkeypatch.context() as context:
            context.setattr(registry, "get_audited_many", altered_metadata)
            context.setattr(study_service, "_make_registry", lambda path: registry)
            with pytest.raises(ValueError, match="durable contract"):
                execution.collect_clustering_run(primary["run_dir"], plan=plan, binding=primary)
