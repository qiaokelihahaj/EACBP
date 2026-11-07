"""Descriptive marker discovery is input-driven and independently audited."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor import ScientificAuditor
from eacbp.auditor.cluster_markers import ClusterMarkersValidator
from eacbp.capabilities import ClusterMarkersCapability, create_default_capability_registry
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import EvidenceType
from eacbp.schemas.task import TaskContract, TaskStatus


def _data(as_sparse=False, n_genes=25):
    counts = np.zeros((8, n_genes), dtype=np.float64)
    counts[:4, 0] = [12, 8, 6, 0]
    counts[4:, 1] = [10, 8, 6, 4]
    counts[:, 2:] = 1
    matrix = sparse.csr_matrix(counts) if as_sparse else counts
    # The clustering representation must never be mistaken for expression.
    return SCData(
        X=matrix * -2, layers={"counts": matrix},
        obs=pd.DataFrame({"cluster": ["A"] * 4 + ["B"] * 4}),
        var=pd.DataFrame({"gene_name": [f"Gene{i:02d}" for i in range(n_genes)]}),
        uns={"normalization": {"log1p": True, "target_sum": 10000.0}},
    )


def _run(tmp_path, data=None, parameters=None):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    in_uri = "adata://markers/annotated/v4"
    registry.register(in_uri, data or _data(), ArtifactType.ANNDATA, "markers", "cluster", "clustering")
    capabilities = create_default_capability_registry()
    contract = TaskContract(
        task_id="task_005b_cluster_markers", capability="cluster_markers", input_artifacts=[in_uri],
        expected_outputs=["table://markers/cluster_markers/v1"], parameters=parameters or {},
        allowed_operations=["summarize_cluster_expression", "rank_cluster_markers"],
        forbidden_operations=["filter_cells", "recluster", "replace_cell_type_labels"],
    )
    capabilities.prepare_contract(contract, registry)
    before = registry.get_metadata(in_uri).sha256_hash
    result = capabilities.execute_contract(contract, registry)
    assert result.status == TaskStatus.SUCCESS, result
    assert registry.get_metadata(in_uri).sha256_hash == before
    return registry, capabilities, contract, result


@pytest.mark.parametrize("as_sparse", [False, True])
def test_markers_rank_existing_clusters_without_a_reference(tmp_path, as_sparse):
    registry, capabilities, contract, result = _run(tmp_path, _data(as_sparse))
    table = registry.get(result.output_artifacts[0])[1]
    assert len(table) == 40
    assert set(table.groupby("cluster").size()) == {20}
    assert table.loc[table["rank"] == 1, "gene"].tolist() == ["Gene00", "Gene01"]
    assert set(table["expression_source"]) == {"layers.counts"}
    assert "p_value" not in table and "fdr_q_value" not in table and "log2fc" not in table
    assert result.inference_contract is None
    report = capabilities.run_validator(contract, result, registry)
    assert report.overall_passed, report
    # Direct suite use must also execute the descriptor's required audit.
    assert ScientificAuditor().audit_task(contract, result, registry).overall_passed
    evidence = capabilities.extract_evidence(contract, result, report, registry)[0]
    assert evidence.type == EvidenceType.CLUSTER_MARKERS
    assert evidence.validation_scope == "within_dataset_descriptive"
    assert evidence.metrics["confirmatory_allowed"] is False
    assert len(evidence.metrics["top_markers_by_cluster"]["A"]) == 20
    assert evidence.inference_contract_id is None


def test_reported_effect_is_mean_log1p_difference_not_fold_change(tmp_path):
    data = _data()
    registry, _, _, result = _run(tmp_path, data)
    table = registry.get(result.output_artifacts[0])[1]
    counts = data.layers["counts"]
    log_expression = np.log1p(counts * (10000.0 / counts.sum(axis=1))[:, None])
    row = table.loc[(table["cluster"] == "A") & (table["gene"] == "Gene00")].iloc[0]
    expected = log_expression[:4, 0].mean() - log_expression[4:, 0].mean()
    assert row["mean_log1p_difference"] == pytest.approx(expected)
    assert row["detection_in"] == 0.75 and row["detection_rest"] == 0
    assert row["effect_units"] == "natural_log1p_normalized_expression"


def test_sparse_markers_and_auditor_never_densify_the_expression_matrix(tmp_path, monkeypatch):
    def prohibit(*args, **kwargs):
        raise AssertionError("Full sparse expression must remain sparse")

    monkeypatch.setattr(sparse.csr_matrix, "toarray", prohibit)
    registry, _, contract, result = _run(tmp_path, _data(True))
    assert ClusterMarkersValidator().audit(contract, result, registry).overall_passed


def test_full_ranking_preserves_duplicate_gene_symbols_and_zero_count_cells(tmp_path):
    data = _data(True)
    data.layers["counts"][3, :] = 0
    data.var.loc[2:3, "gene_name"] = "DuplicatedSymbol"
    registry, _, contract, result = _run(tmp_path, data, {"top_n": None})
    table = registry.get(result.output_artifacts[0])[1]
    assert len(table) == 2 * data.n_vars
    assert not table.duplicated(["cluster", "gene_index"]).any()
    assert len(table[table["gene"] == "DuplicatedSymbol"]) == 4
    assert ClusterMarkersValidator().audit(contract, result, registry).overall_passed


@pytest.mark.parametrize("normalized", [False, True])
def test_explicit_normalized_x_or_integer_counts_are_supported(tmp_path, normalized):
    counts = np.array([[4, 0], [2, 0], [0, 3], [0, 2]], dtype=np.float64)
    expression = np.log1p(counts * (10000 / counts.sum(axis=1))[:, None]) if normalized else counts
    data = SCData(expression, pd.DataFrame({"cluster": ["A", "A", "B", "B"]}),
                  pd.DataFrame(index=["alpha", "beta"]),
                  uns={"normalization": {"log1p": True, "target_sum": 10000.0}} if normalized else {})
    registry, _, contract, result = _run(tmp_path, data)
    table = registry.get(result.output_artifacts[0])[1]
    assert len(table) == 4
    assert set(table["expression_source"]) == {"X_normalized" if normalized else "X_integer_counts"}
    assert ClusterMarkersValidator().audit(contract, result, registry).overall_passed


@pytest.mark.parametrize("corruption", ["detection", "gene", "rank", "missing_cluster", "missing_gene", "p_value", "effect", "source", "metrics"])
def test_independent_audit_rejects_incorrect_values_ranks_coverage_and_inference(tmp_path, corruption):
    registry, _, contract, result = _run(tmp_path)
    meta, original = registry.get(result.output_artifacts[0])
    table = original.copy()
    if corruption == "detection":
        table.loc[0, "detection_in"] = 0.2
    elif corruption == "gene":
        table.loc[0, "gene"] = "FabricatedGene"
    elif corruption == "rank":
        table.loc[0, "gene_index"] = 24
    elif corruption == "missing_cluster":
        table = table[table["cluster"] == "A"]
    elif corruption == "missing_gene":
        table = table.iloc[1:]
    elif corruption == "p_value":
        table["p_value"] = 0.001
    elif corruption == "effect":
        table.loc[0, "mean_log1p_difference"] += 1
    elif corruption == "source":
        table["expression_source"] = "X_normalized"
    elif corruption == "metrics":
        result.metrics["confirmatory_allowed"] = True
    uri = "table://markers/corrupted_markers/v1"
    registry.register(uri, table, ArtifactType.TABLE, "markers", contract.task_id, meta.operation,
                      parent_uris=meta.parent_uris, parameters=meta.parameters, summary_metrics=meta.summary_metrics)
    corrupted = result.model_copy(update={"output_artifacts": [uri]})
    report = ClusterMarkersValidator().audit(contract, corrupted, registry)
    assert not report.overall_passed, report
    assert any(check.check_name == "cluster_marker_integrity" and not check.passed for check in report.checks)


@pytest.mark.parametrize("problem", ["single_cluster", "missing_cluster", "nan", "negative", "undocumented_transform"])
def test_marker_discovery_rejects_invalid_or_ambiguous_inputs(tmp_path, problem):
    data = _data()
    if problem == "single_cluster":
        data.obs["cluster"] = "A"
    elif problem == "missing_cluster":
        data.obs.loc[0, "cluster"] = None
    elif problem == "nan":
        data.layers["counts"][0, 0] = np.nan
    elif problem == "negative":
        data.layers["counts"][0, 0] = -1
    elif problem == "undocumented_transform":
        data.X = np.abs(data.X) + 0.25
        data.layers.clear()
        data.uns.clear()
    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://markers/annotated/v4"
    registry.register(uri, data, ArtifactType.ANNDATA, "markers", "cluster", "clustering")
    contract = TaskContract(task_id="markers", capability="cluster_markers", input_artifacts=[uri])
    with pytest.raises(ValueError):
        ClusterMarkersCapability().execute(contract, registry)


def test_descriptor_contract_scope_and_plan_use_one_shared_task():
    capability = ClusterMarkersCapability()
    descriptor = capability.get_descriptor()
    assert descriptor.scope == "shared"
    assert descriptor.required_audit_ids == ("cluster_marker_integrity",)
    with pytest.raises(ValidationError):
        descriptor.validate_parameters({"top_n": 19})
    cluster = TaskContract(task_id="cluster", capability="clustering", expected_outputs=["adata://s/annotated/v4"])
    task = descriptor.build_plan(manifest=SimpleNamespace(study_id="s"), parameters={}, tasks=[cluster])
    assert task.input_artifacts == ["adata://s/annotated/v4"]
    assert task.expected_outputs == ["table://s/cluster_markers/v1"]
    assert task.depends_on == ["cluster"]
    assert descriptor.build_plan(manifest=SimpleNamespace(study_id="s"), parameters={}, tasks=[cluster, task]) == []
