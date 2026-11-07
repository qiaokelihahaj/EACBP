"""Native velocity API contracts, count preservation, and independent audit."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor import ScientificAuditor
from eacbp.auditor.rna_velocity import RNAVelocityValidator
from eacbp.capabilities import RNAVelocityCapability, create_default_capability_registry
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import EvidenceType
from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest
from eacbp.schemas.task import ExecutionFailureType, TaskContract, TaskStatus


def _input():
    cells = [f"cell{i}" for i in range(12)]
    spliced = sparse.csr_matrix(np.arange(48).reshape(12, 4) % 8 + 1, dtype=np.float32)
    unspliced = sparse.csr_matrix(np.arange(48).reshape(12, 4) % 3 + 1, dtype=np.float32)
    edges = sparse.csr_matrix(np.ones((12, 12)) - np.eye(12), dtype=np.float32)
    return SCData(
        X=spliced * 0.25, obs=pd.DataFrame({"cell_id": cells, "cluster": ["0"] * 6 + ["1"] * 6}, index=cells),
        var=pd.DataFrame({"gene_name": ["a", "b", "c", "d"]}, index=["a", "b", "c", "d"]),
        layers={"counts": spliced, "spliced": spliced, "unspliced": unspliced},
        obsm={"X_pca": np.arange(36, dtype=np.float32).reshape(12, 3)},
        obsp={"connectivities": edges, "distances": edges * 2},
        uns={"clustering": {"method": "scanpy_leiden_umap_v1"},
             "neighbors": {"connectivities_key": "connectivities", "distances_key": "distances", "params": {"n_neighbors": 12}}},
    )


class FakeNativeScVelo:
    """API fixture, not a scientific reproduction or substitute method."""
    __version__ = "0.3.4"

    def __init__(self):
        self.calls = []
        self.pp = SimpleNamespace(normalize_per_cell=self.normalize, moments=self.moments)
        self.tl = SimpleNamespace(velocity=self.velocity, velocity_graph=self.graph,
                                  terminal_states=self.terminal, velocity_pseudotime=self.pseudotime,
                                  velocity_confidence=self.confidence)

    def normalize(self, adata, **kwargs):
        self.calls.append(("normalize", kwargs))
        assert kwargs["enforce"] is True and kwargs["use_initial_size"] is False
        assert kwargs["layers"] == ["spliced", "unspliced"]
        assert (adata.X - adata.layers["counts"]).nnz == 0
        for layer in ["X", "spliced", "unspliced"]:
            raw = adata.X if layer == "X" else adata.layers[layer]
            totals = np.asarray(raw.sum(axis=1)).ravel()
            normalized = raw.multiply((kwargs["counts_per_cell_after"] / np.where(totals > 0, totals, 1))[:, None]).tocsr()
            if layer == "X":
                adata.X = normalized
            else:
                adata.layers[layer] = normalized

    def moments(self, adata, **kwargs):
        self.calls.append(("moments", kwargs))
        assert kwargs["n_neighbors"] is None and kwargs["use_highly_variable"] is False
        graph = (adata.obsp["connectivities"] > 0).astype(np.float32)
        graph.setdiag(1)
        graph = graph.multiply(1 / np.asarray(graph.sum(axis=1))).tocsr().astype(np.float32)
        adata.layers["Ms"] = (graph @ adata.layers["spliced"]).toarray().astype(np.float32)
        adata.layers["Mu"] = (graph @ adata.layers["unspliced"]).toarray().astype(np.float32)

    def velocity(self, adata, **kwargs):
        self.calls.append(("velocity", kwargs))
        assert kwargs == {"mode": "stochastic", "filter_genes": False, "use_highly_variable": None}
        adata.layers["velocity"] = adata.layers["Mu"] - adata.layers["Ms"] * 0.5
        adata.var["velocity_genes"] = True
        adata.var["velocity_gamma"] = 0.5
        adata.var["velocity_r2"] = 0.6
        adata.uns["velocity_params"] = {"mode": "stochastic", "fit_offset": False}

    def graph(self, adata, **kwargs):
        self.calls.append(("graph", kwargs))
        assert kwargs["n_jobs"] == 1
        adata.uns["velocity_graph"] = sparse.diags(np.ones(adata.n_obs - 1) * 0.8, 1, shape=(adata.n_obs, adata.n_obs), format="csr")
        adata.uns["velocity_graph_neg"] = sparse.diags(np.ones(adata.n_obs - 1) * -0.2, -1, shape=(adata.n_obs, adata.n_obs), format="csr")

    def terminal(self, adata, **kwargs):
        self.calls.append(("terminal", kwargs))
        assert kwargs["groupby"] is None
        assert "cell_fate" not in adata.obs and "iroot" not in adata.uns
        assert "root_cells" not in adata.obs and "end_points" not in adata.obs
        adata.obs["root_cells"] = np.linspace(1, 0, adata.n_obs)
        adata.obs["end_points"] = np.linspace(0, 1, adata.n_obs)

    def pseudotime(self, adata, **kwargs):
        self.calls.append(("pseudotime", kwargs))
        assert kwargs["groupby"] is None
        assert kwargs["root_key"] == "root_cells" and kwargs["end_key"] == "end_points"
        adata.obs["velocity_pseudotime"] = np.linspace(0, 1, adata.n_obs)

    def confidence(self, adata):
        self.calls.append(("confidence", {}))
        adata.obs["velocity_confidence"] = 0.6


def _mock_native(monkeypatch, native):
    import eacbp.capabilities.rna_velocity as module
    real_import = module.importlib.import_module

    def load(name, *args, **kwargs):
        if name == "scvelo":
            return native
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(module.importlib, "import_module", load)
    monkeypatch.setattr(module.metadata, "version", lambda package: "0.3.4")
    # This API fixture does no scVelo generalized least-squares arithmetic.
    # Actual library compatibility is covered by the real integration gate.
    monkeypatch.setattr(module, "_require_velocity_numpy", lambda: None)


def _run(tmp_path, monkeypatch, data=None, native=None):
    pytest.importorskip("anndata")
    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://velocity/annotated/v4"
    registry.register(uri, data or _input(), ArtifactType.ANNDATA, "velocity", "cluster", "scanpy_leiden_umap_v1")
    native = native or FakeNativeScVelo()
    _mock_native(monkeypatch, native)
    capabilities = create_default_capability_registry()
    contract = TaskContract(task_id="task_005c_rna_velocity", capability="rna_velocity", input_artifacts=[uri],
                            expected_outputs=["adata://velocity/rna_velocity/v1", "table://velocity/rna_velocity_summary/v1"],
                            allowed_operations=list(RNAVelocityCapability.contract_operations),
                            forbidden_operations=["filter_cells", "filter_genes", "recluster"])
    result = capabilities.execute_contract(contract, registry)
    assert result.status == TaskStatus.SUCCESS, result
    return registry, capabilities, contract, result, native


def test_native_api_pipeline_and_independent_summary_audit(tmp_path, monkeypatch):
    registry, capabilities, contract, result, native = _run(tmp_path, monkeypatch)
    assert [name for name, _ in native.calls] == ["normalize", "moments", "velocity", "graph", "terminal", "pseudotime", "confidence"]
    output = registry.get_scdata(result.output_artifacts[0])[1]
    assert output.shape == (12, 4)
    assert result.inference_contract is None and result.scientific_result is None
    report = capabilities.run_validator(contract, result, registry)
    assert report.overall_passed, report
    assert ScientificAuditor().audit_task(contract, result, registry).overall_passed
    evidence = capabilities.extract_evidence(contract, result, report, registry)[0]
    assert evidence.type == EvidenceType.RNA_VELOCITY
    assert evidence.validation_scope == "within_dataset_model_inferred"
    assert evidence.metrics["confirmatory_allowed"] is False
    assert output.uns["rna_velocity"]["kinetic_groupby"] == "none_global_graph"
    assert output.uns["rna_velocity"]["caller_root_used"] is False


def test_existing_model_roots_and_fates_are_cleared_before_inference(tmp_path, monkeypatch):
    data = _input()
    data.obs["cell_fate"] = "published answer"
    data.obs["root_cells"] = 1.0
    data.obs["end_points"] = 0.0
    data.obs["latent_time"] = 100.0
    data.uns["iroot"] = 10
    data.var["spearmans_score"] = 0
    registry, _, contract, result, _ = _run(tmp_path, monkeypatch, data)
    output = registry.get_scdata(result.output_artifacts[0])[1]
    assert "cell_fate" not in output.obs and "latent_time" not in output.obs and "iroot" not in output.uns
    assert output.obs["root_cells"].iloc[0] == 1 and output.obs["root_cells"].iloc[-1] == 0
    assert RNAVelocityValidator().audit(contract, result, registry).overall_passed


def test_constant_gamma_omission_is_an_explicit_audit_limit(tmp_path, monkeypatch):
    class ConstantGammaNative(FakeNativeScVelo):
        def velocity(self, adata, **kwargs):
            super().velocity(adata, **kwargs)
            del adata.var["velocity_gamma"]

    registry, _, contract, result, _ = _run(tmp_path, monkeypatch, native=ConstantGammaNative())
    report = RNAVelocityValidator().audit(contract, result, registry)
    assert report.overall_passed, report
    check = next(check for check in report.checks if check.check_name == "rna_velocity_parameter_and_residual_consistency")
    assert check.metrics["residual_equation_evaluated"] is False
    assert check.metrics["kinetic_fit_refitted"] is False


@pytest.mark.parametrize("corruption", ["counts", "raw_spliced", "normalized", "moments", "velocity", "finite_velocity", "gamma", "r2",
                                     "graph_shape", "graph_range", "summary", "missing_cell", "caller_root", "dynamical", "cluster", "pca",
                                     "parameters", "seed", "n_dcs", "target_sum", "kinetic_groupby"])
def test_independent_audit_rejects_corrupted_kinetics_counts_graph_and_summary(tmp_path, monkeypatch, corruption):
    registry, _, contract, result, _ = _run(tmp_path, monkeypatch)
    meta, payload = registry.get(result.output_artifacts[0])
    data = SCData.from_payload(payload).copy()
    table_meta, table = registry.get(result.output_artifacts[1])
    table = table.copy()
    if corruption == "counts":
        data.layers["counts"][0, 0] += 1
    elif corruption == "raw_spliced":
        data.layers["spliced_counts"][0, 0] += 1
    elif corruption == "normalized":
        data.layers["spliced"][0, 0] += 100
    elif corruption == "moments":
        data.layers["Ms"][0, 0] += 100
    elif corruption == "velocity":
        data.layers["velocity"][0, 0] = np.nan
    elif corruption == "finite_velocity":
        data.layers["velocity"][0, 0] += 100
    elif corruption == "gamma":
        data.var.loc["a", "velocity_gamma"] = 0.8
    elif corruption == "r2":
        data.var.loc["a", "velocity_r2"] = np.nan
    elif corruption == "graph_shape":
        data.uns["velocity_graph"] = sparse.eye(3, format="csr")
    elif corruption == "graph_range":
        data.uns["velocity_graph"][0, 1] = 2.0
    elif corruption == "summary":
        table.loc[0, "root_score"] = 0.25
    elif corruption == "missing_cell":
        table = table.iloc[1:]
    elif corruption == "caller_root":
        data.uns["rna_velocity"]["caller_root_used"] = True
    elif corruption == "dynamical":
        data.uns["rna_velocity"]["mode"] = "dynamical"
    elif corruption == "cluster":
        data.obs.loc["cell0", "cluster"] = "1"
    elif corruption == "pca":
        data.obsm["X_pca"][0, 0] += 1
    elif corruption == "parameters":
        meta.parameters["target_sum"] = 7777
    elif corruption == "n_dcs":
        data.uns["rna_velocity"]["n_dcs"] = 3
    elif corruption == "target_sum":
        data.uns["rna_velocity"]["target_sum"] = 7777
    elif corruption == "kinetic_groupby":
        data.uns["rna_velocity"]["kinetic_groupby"] = "cluster"
    out_uri, table_uri = "adata://velocity/corrupt_velocity/v1", "table://velocity/corrupt_summary/v1"
    registry.register(out_uri, data, ArtifactType.ANNDATA, "velocity", contract.task_id, meta.operation,
                      parent_uris=meta.parent_uris, parameters=meta.parameters, software_versions=meta.software_versions,
                      summary_metrics=meta.summary_metrics, random_seed=999 if corruption == "seed" else meta.random_seed)
    registry.register(table_uri, table, ArtifactType.TABLE, "velocity", contract.task_id, table_meta.operation,
                      parent_uris=[out_uri], parameters=table_meta.parameters, software_versions=table_meta.software_versions,
                      summary_metrics=table_meta.summary_metrics, random_seed=table_meta.random_seed)
    corrupted = result.model_copy(update={"output_artifacts": [out_uri, table_uri]})
    report = RNAVelocityValidator().audit(contract, corrupted, registry)
    assert not report.overall_passed, report
    assert any(check.check_name == "rna_velocity_integrity" and not check.passed for check in report.checks)


@pytest.mark.parametrize("problem", ["missing_layer", "fractional_counts", "negative_counts", "missing_neighbors", "non_native_cluster"])
def test_native_velocity_rejects_unqualified_inputs_before_library_import(tmp_path, monkeypatch, problem):
    data = _input()
    if problem == "missing_layer":
        del data.layers["unspliced"]
    elif problem == "fractional_counts":
        data.layers["unspliced"][0, 0] = 0.5
    elif problem == "negative_counts":
        data.layers["spliced"][0, 0] = -1
    elif problem == "missing_neighbors":
        data.obsp.clear()
    elif problem == "non_native_cluster":
        data.uns["clustering"]["method"] = "kmeans_marker_embedding_v1"
    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://velocity/annotated/v4"
    registry.register(uri, data, ArtifactType.ANNDATA, "velocity", "cluster", "clustering")
    with pytest.raises(ValueError):
        RNAVelocityCapability().execute(TaskContract(task_id="v", capability="rna_velocity", input_artifacts=[uri]), registry)


def test_missing_scvelo_is_a_dependency_failure_without_fallback(tmp_path, monkeypatch):
    import eacbp.capabilities.rna_velocity as module
    real_import = module.importlib.import_module

    def missing(name, *args, **kwargs):
        if name == "scvelo":
            raise ImportError("missing fixture dependency")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(module.importlib, "import_module", missing)
    monkeypatch.setattr(module, "_require_velocity_numpy", lambda: None)
    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://velocity/annotated/v4"
    registry.register(uri, _input(), ArtifactType.ANNDATA, "velocity", "cluster", "clustering")
    result = create_default_capability_registry().execute_contract(TaskContract(task_id="v", capability="rna_velocity", input_artifacts=[uri]), registry)
    assert result.status == TaskStatus.EXECUTION_FAILURE
    assert result.error_type == ExecutionFailureType.DEPENDENCY_ERROR
    assert "no substitute" in result.error_message
    assert not result.output_artifacts


@pytest.mark.parametrize("numpy_version", ["2.4.0", "2.5.3"])
def test_unsupported_numpy_is_an_actionable_dependency_error(numpy_version, monkeypatch):
    import eacbp.capabilities.rna_velocity as module
    monkeypatch.setattr(module.np, "__version__", numpy_version)
    with pytest.raises(ImportError, match=r"requires numpy>=1.24,<2.4"):
        module._require_velocity_numpy()


def test_supported_numpy_version_passes_compatibility_check(monkeypatch):
    import eacbp.capabilities.rna_velocity as module
    monkeypatch.setattr(module.np, "__version__", "2.3.5")
    module._require_velocity_numpy()


@pytest.mark.parametrize("parameters", [{"root_cell_id": "known"}, {"root_key": "published_root"}, {"mode": "dynamical"}, {"groupby": "clusters"}, {"n_jobs": 2}])
def test_descriptor_rejects_external_orientation_and_parallel_override(parameters):
    with pytest.raises(ValidationError):
        RNAVelocityCapability().get_descriptor().validate_parameters(parameters)


def test_descriptor_plans_once_after_native_clustering():
    descriptor = RNAVelocityCapability().get_descriptor()
    assert descriptor.scope == "shared"
    task = TaskContract(task_id="c", capability="clustering", method="scanpy_leiden_umap_v1", expected_outputs=["adata://s/annotated/v4"])
    velocity = descriptor.build_plan(manifest=SimpleNamespace(study_id="s"), parameters={}, tasks=[task])
    assert velocity.task_id == "task_005c_rna_velocity"
    assert velocity.depends_on == ["c"]
    assert descriptor.build_plan(manifest=SimpleNamespace(study_id="s"), parameters={}, tasks=[task, velocity]) == []
    task.method = "kmeans_marker_embedding_v1"
    with pytest.raises(ValueError, match="native"):
        descriptor.build_plan(manifest=SimpleNamespace(study_id="s"), parameters={}, tasks=[task])
    task.method = None
    assert descriptor.build_plan(manifest=SimpleNamespace(study_id="s"), parameters={}, tasks=[task]).capability == "rna_velocity"


def test_native_autonomous_dag_build_resolves_extension_context_without_manual_clustering_method():
    from eacbp.orchestrator.dag import ComputationalDAGPlanner
    from eacbp.orchestrator.planning import build_study_tasks, resolve_task_contract

    registry = create_default_capability_registry()
    manifest = StudyManifest(study_id="native_velocity_dag", biological_design=BiologicalDesign(species="mouse", tissue="test"),
                             data=DataSpec(raw_artifact_uri="adata://native_velocity_dag/raw/v1", has_rna_velocity=True))
    state = {"method_profile": "standard", "autonomous_analysis": {"intent": "exploration"},
             "analysis_extensions": {"rna_velocity": True, "cluster_markers": True}}
    tasks = build_study_tasks(manifest, state, registry)
    abstract_tasks = ComputationalDAGPlanner.build_study_plan(
        manifest, {key: value for key, value in state.items() if key != "analysis_extensions"}, expand_targets=False)
    original_contracts = {task.task_id: task.model_dump() for task in abstract_tasks}
    assert {task.task_id: task.model_dump() for task in tasks if task.task_id in original_contracts} == original_contracts
    clustering = next(task for task in tasks if task.capability == "clustering")
    selected = [task for task in tasks if task.capability == "rna_velocity"]
    assert len(selected) == 1
    velocity = selected[0]
    assert velocity.depends_on == [clustering.task_id]
    assert velocity.input_artifacts == clustering.expected_outputs
    assert velocity.expected_outputs == ["adata://native_velocity_dag/rna_velocity/v1", "table://native_velocity_dag/rna_velocity_summary/v1"]
    resolve_task_contract(clustering, manifest, state, registry)
    assert clustering.method == "scanpy_leiden_umap_v1"
    resolve_task_contract(velocity, manifest, state, registry)
    assert velocity.method == "scvelo_stochastic_v1"
    assert "rna_velocity_integrity" in velocity.validation_requirements
    assert "target_branch" not in velocity.parameters
    with pytest.raises(ValueError, match="native"):
        build_study_tasks(manifest, {**state, "method_overrides": {"clustering": "kmeans_marker_embedding_v1"}}, registry)
    with pytest.raises(ValueError, match="native"):
        build_study_tasks(manifest, {**state, "method_profile": "baseline"}, registry)


@pytest.mark.parametrize("option", ["fit_offset", "min_r2", "perc", "use_raw", "groups", "gene_subset"])
def test_unsupported_scientific_options_cannot_be_recorded_as_executed(option):
    with pytest.raises(ValidationError, match="unimplemented scientific overrides"):
        RNAVelocityCapability().get_descriptor().validate_parameters({option: True})
