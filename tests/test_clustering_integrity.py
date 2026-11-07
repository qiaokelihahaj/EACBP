"""Real Leiden is reconstructed independently after durable serialization."""

from copy import deepcopy
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

pytest.importorskip("scanpy")
pytest.importorskip("leidenalg")

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor import ScientificAuditor
from eacbp.auditor.clustering import ClusteringIntegrityValidator, validate_clustering_execution_receipt
from eacbp.capabilities.integration import IntegrationCapability, compute_pca
from eacbp.capabilities.real_methods import LeidenClusteringCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.audit import model_to_jsonable, stable_fingerprint
from eacbp.schemas.task import TaskContract


def _source(as_sparse=True, counts_layer=True, n_features=18, n_hvg=14):
    rng = np.random.default_rng(31)
    counts = rng.poisson(1.2, (72, n_features)).astype(np.int64)
    for group in range(3):
        counts[group * 24:(group + 1) * 24, group * 4:(group + 1) * 4] += 8
    x = sparse.csr_matrix(np.log1p(counts).astype(np.float32)) if as_sparse else np.log1p(counts).astype(np.float32)
    layer = sparse.csr_matrix(counts) if as_sparse else counts
    return SCData(x,
                  pd.DataFrame({"cell_id": [f"c{i}" for i in range(72)], "donor_id": [f"d{i % 3}" for i in range(72)],
                                "batch_id": [f"b{i % 2}" for i in range(72)]}, index=[f"c{i}" for i in range(72)]),
                  pd.DataFrame({"gene_name": [f"g{i}" for i in range(n_features)],
                                "highly_variable": [True] * n_hvg + [False] * (n_features - n_hvg)},
                               index=[f"g{i}" for i in range(n_features)]),
                  layers={"counts": layer} if counts_layer else {})


def _execute(directory, *, integrated=True, as_sparse=True, counts_layer=True):
    registry = ArtifactRegistry(directory)
    uri = "adata://cluster/normalized/v2"
    data = _source(as_sparse=as_sparse, counts_layer=counts_layer)
    registry.register(uri, data, ArtifactType.ANNDATA, "cluster", "normalization", "normalization")
    if integrated:
        task = TaskContract(task_id="integration", capability="integration", method="no_correction_v1",
                            input_artifacts=[uri], parameters={"n_components": 5, "random_seed": 17})
        integrated_result = IntegrationCapability("no_correction_v1").execute(task, registry)
        uri = integrated_result.output_artifacts[0]
    contract = TaskContract(task_id="clustering", capability="clustering", method="scanpy_leiden_umap_v1",
                            input_artifacts=[uri], expected_outputs=["adata://cluster/annotated/v4"],
                            parameters={"n_components": 7, "n_neighbors": 10, "resolution": 1.0, "random_seed": 13},
                            validation_requirements=["leiden_execution_integrity"])
    result = LeidenClusteringCapability().execute(contract, registry)
    return registry, contract, result


@pytest.fixture(scope="module")
def native(tmp_path_factory):
    return _execute(tmp_path_factory.mktemp("native_clustering"))


def test_reopened_native_has_effective_settings_and_independent_audit(native):
    registry, contract, result = native
    reopened = ArtifactRegistry(registry.storage.base_dir)
    report = ClusteringIntegrityValidator().audit(contract, result, reopened)
    assert report.overall_passed, report
    receipt = validate_clustering_execution_receipt(contract, result, reopened)
    assert receipt["representation"]["reused"] is True
    assert receipt["representation"]["requested_n_components"] == 7
    assert receipt["representation"]["actual_n_components"] == 5
    assert receipt["representation"]["algorithm"] == "uncentered_sparse_truncated_svd"
    assert receipt["representation"]["feature_selection"] == "highly_variable"
    assert receipt["representation"]["inherited_representation_receipt"]["n_features_used"] == 14
    assert receipt["neighbors"]["all_representation_components_used"] is True
    assert receipt["input_roles"]["donor_independence_established"] is False
    assert receipt["independent_biological_validation"] is False
    assert report.checks[0].metrics["partition_recomputed"] is True
    assert report.checks[0].metrics["umap_recomputed"] is False
    # Direct suite callers may omit descriptor preparation but cannot omit the independent audit.
    direct = contract.model_copy(update={"validation_requirements": []})
    suite = ScientificAuditor().audit_task(direct, result, reopened)
    assert suite.overall_passed, suite
    assert any(c.check_name == "leiden_execution_integrity" and c.passed for c in suite.checks)


@pytest.mark.parametrize("as_sparse,counts_layer", [(False, True), (True, False)])
def test_missing_input_pca_is_reconstructed_using_all_features(tmp_path, as_sparse, counts_layer):
    registry, contract, result = _execute(tmp_path / "missing_pca", integrated=False,
                                          as_sparse=as_sparse, counts_layer=counts_layer)
    report = ClusteringIntegrityValidator().audit(contract, result, registry)
    assert report.overall_passed, report
    receipt = result.metrics["execution_receipt"]
    assert receipt["representation"]["reused"] is False
    assert receipt["representation"]["feature_selection"] == "all_features"
    assert receipt["representation"]["algorithm"] == "scanpy_centered_pca_arpack"
    assert receipt["representation"]["actual_n_components"] == 7
    assert receipt["counts_source"] == ("layers.counts" if counts_layer else "X_input_representation")


@pytest.mark.parametrize("route", ["sparse", "dense", "fallback"])
def test_no_correction_receipt_records_original_numerical_route(tmp_path, monkeypatch, route):
    data = _source(as_sparse=route != "dense")
    x = data.X[:, data.var.highly_variable.to_numpy()]
    seed, components = 17, 5
    if route == "sparse":
        from scipy.sparse.linalg import svds
        v0 = np.random.default_rng(seed).standard_normal(min(x.shape))
        u, singular, _ = svds(x.astype(np.float32), k=components, v0=v0)
        original = u[:, ::-1] * singular[::-1]
    else:
        if route == "fallback":
            import scipy.sparse.linalg
            monkeypatch.setattr(scipy.sparse.linalg, "svds", lambda *a, **k: (_ for _ in ()).throw(ValueError("forced fallback")))
            dense = x.toarray()
        else:
            dense = np.asarray(x, dtype=np.float32)
        centered = dense - np.mean(dense, axis=0, keepdims=True)
        u, singular, _ = np.linalg.svd(centered, full_matrices=False)
        original = u[:, :components] * singular[:components]
    assert np.array_equal(compute_pca(x, components, seed), original)
    registry = ArtifactRegistry(tmp_path / "route")
    uri = "adata://route/normalized/v2"
    registry.register(uri, data, ArtifactType.ANNDATA, "route", "norm", "norm")
    task = TaskContract(task_id="integration", capability="integration", input_artifacts=[uri],
                        parameters={"n_components": components, "random_seed": seed})
    result = IntegrationCapability("no_correction_v1").execute(task, registry)
    output = registry.load_scdata(result.output_artifacts[0])
    assert np.array_equal(output.obsm["X_pca"], original)
    receipt = output.uns["integration"]["representation_receipt"]
    assert receipt["algorithm"] == {"sparse": "uncentered_sparse_truncated_svd", "dense": "centered_dense_svd",
                                    "fallback": "centered_dense_svd_fallback"}[route]
    assert receipt["zero_center"] is (route != "sparse")
    assert receipt["fallback_used"] is (route == "fallback")


def _rehashed(native, tmp_path, corruption):
    original, contract, result = native
    replacement = ArtifactRegistry(tmp_path / "rehashed")
    # Preserve the already verified input/ancestor files; only the attacked output
    # is published afresh, with a valid newly computed payload hash.
    for meta in original.list_artifacts():
        if meta.uri not in result.output_artifacts:
            copied = replacement.storage.base_dir / Path(meta.storage_path).relative_to(original.storage.base_dir)
            copied.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(meta.storage_path, copied)
            replacement.registry[meta.uri] = meta.model_copy(update={"storage_path": str(copied)}, deep=True)
    replacement._rebuild_lineage()
    replacement._persist_index()
    meta, output = original.get_scdata(result.output_artifacts[0])
    output = output.copy()
    parameters = deepcopy(meta.parameters)
    parents = list(meta.parent_uris)
    metrics = deepcopy(meta.summary_metrics)
    changed = result.model_copy(deep=True)
    receipt = output.uns["clustering"]["execution_receipt"]
    if corruption == "counts":
        output.layers["counts"].data[0] += 1
    elif corruption == "expression":
        output.X.data[0] += .1
    elif corruption == "cell_order":
        output = output.subset_obs(np.arange(output.n_obs)[::-1])
    elif corruption == "feature_metadata":
        output.var.iloc[0, output.var.columns.get_loc("gene_name")] = "forged_gene"
    elif corruption == "feature_mask":
        output.var.iloc[-1, output.var.columns.get_loc("highly_variable")] = True
    elif corruption == "donor":
        output.obs["donor_id"] = output.obs["donor_id"].astype(object)
        output.obs.iloc[0, output.obs.columns.get_loc("donor_id")] = "invented_donor"
    elif corruption == "pca":
        output.obsm["X_pca"][0, 0] += .1
    elif corruption in {"connectivities", "distances"}:
        output.obsp[corruption].data[0] += .1
    elif corruption == "neighbors_uns":
        output.uns["neighbors"]["params"]["n_neighbors"] += 1
    elif corruption == "leiden_uns":
        output.uns["cluster"]["params"]["resolution"] = 3.0
    elif corruption in {"partition_coordinated", "label_rename"}:
        labels = output.obs["cluster"].astype(str).to_numpy()
        if corruption == "partition_coordinated":
            labels[0] = next(label for label in sorted(set(labels)) if label != labels[0])
        else:
            labels = np.array(["renamed_" + label for label in labels])
        output.obs["cluster"] = labels
        output.obs["leiden"] = labels.copy()
        output.obs["cell_type"] = [f"Cluster_{label}" for label in labels]
        codes, _ = pd.factorize(labels, sort=False)
        receipt["leiden"]["partition_sha256"] = stable_fingerprint({"cell_ids": output.obs.index.astype(str).tolist(), "partition": codes.tolist()})
        metrics = deepcopy(model_to_jsonable(output.uns["clustering"]))
        metrics.pop("method")
        changed.metrics = deepcopy(metrics)
    elif corruption == "label_alias":
        output.obs["leiden"] = "wrong"
    elif corruption == "cell_type":
        output.obs["cell_type"] = "fabricated_identity"
    elif corruption == "source_sha":
        receipt["source_sha256"] = "sha256:" + "0" * 64
    elif corruption == "receipt_extra":
        receipt["p_value"] = .0001
    elif corruption == "receipt_resolution":
        receipt["leiden"]["resolution"] = 2.0
    elif corruption == "receipt_components":
        receipt["representation"]["actual_n_components"] = 7
    elif corruption == "receipt_algorithm":
        receipt["representation"]["algorithm"] = "centered_dense_svd"
    elif corruption == "receipt_partition":
        receipt["leiden"]["partition_sha256"] = "0" * 64
    elif corruption == "receipt_inherited_mask":
        receipt["representation"]["inherited_representation_receipt"]["feature_mask_sha256"] = "0" * 64
    elif corruption == "receipt_independence":
        receipt["input_roles"]["donor_independence_established"] = True
    elif corruption == "metadata_resolution":
        parameters["resolution"] = 2.0
    elif corruption == "metadata_parents":
        parents = ["adata://cluster/normalized/v2"]
    elif corruption == "metadata_summary":
        metrics["n_clusters"] += 1
    elif corruption == "result_receipt":
        changed.metrics["execution_receipt"]["leiden"]["resolution"] = 2.0
    elif corruption == "unknown_option":
        parameters["n_pcs"] = 3
        contract = contract.model_copy(update={"parameters": parameters})
    else:
        raise AssertionError(corruption)
    replacement.register(meta.uri, output, meta.type, meta.study_id, meta.created_by_task, meta.operation,
                         parent_uris=parents, parameters=parameters, summary_metrics=metrics, random_seed=meta.random_seed,
                         software_versions=meta.software_versions)
    replacement.verify_many(changed.output_artifacts)
    return replacement, contract, changed


@pytest.mark.parametrize("corruption", [
    "counts", "expression", "cell_order", "feature_metadata", "feature_mask", "donor", "pca",
    "connectivities", "distances", "neighbors_uns", "leiden_uns", "partition_coordinated", "label_alias", "cell_type",
    "source_sha", "receipt_extra", "receipt_resolution", "receipt_components", "receipt_algorithm", "receipt_partition",
    "receipt_inherited_mask", "receipt_independence", "metadata_resolution", "metadata_parents", "metadata_summary",
    "result_receipt", "unknown_option",
])
def test_valid_rehash_does_not_admit_corrupted_execution(native, tmp_path, corruption):
    registry, contract, result = _rehashed(native, tmp_path, corruption)
    report = ClusteringIntegrityValidator().audit(contract, result, registry)
    assert not report.overall_passed, (corruption, report)
    assert report.checks[0].check_name == "leiden_execution_integrity"


def test_partition_equivalence_allows_cluster_name_permutation(native, tmp_path):
    registry, contract, result = _rehashed(native, tmp_path, "label_rename")
    report = ClusteringIntegrityValidator().audit(contract, result, registry)
    assert report.overall_passed, report


def test_cluster_identifiers_must_not_collide_after_conversion_to_strings(native, monkeypatch):
    registry, contract, result = native
    original = registry.get_scdata
    def changed(uri):
        meta, data = original(uri)
        if uri == result.output_artifacts[0]:
            labels = data.obs["cluster"].astype(object)
            labels.iloc[0] = int(labels.iloc[0])
            data.obs["cluster"] = labels
            data.obs["leiden"] = labels.copy()
        return meta, data
    monkeypatch.setattr(registry, "get_scdata", changed)
    with pytest.raises(ValueError, match="Cluster labels"):
        validate_clustering_execution_receipt(contract, result, registry)


def test_presenting_admitted_receipt_does_not_refit(native, monkeypatch):
    import scanpy as sc
    registry, contract, result = native
    def forbidden(*args, **kwargs):
        raise AssertionError("Receipt presentation must not refit")
    monkeypatch.setattr(sc.pp, "neighbors", forbidden)
    monkeypatch.setattr(sc.pp, "pca", forbidden)
    monkeypatch.setattr(sc.tl, "leiden", forbidden)
    receipt = validate_clustering_execution_receipt(contract, result, registry, recompute_partition=False)
    assert receipt["validation_scope"] == "execution_reconstruction_only"
    assert receipt["umap_reconstruction_in_scope"] is False


def test_native_task_executor_transaction_and_durable_admission(tmp_path):
    from eacbp.capabilities import create_default_capability_registry
    from eacbp.orchestrator.execution import TaskExecutor

    registry = ArtifactRegistry(tmp_path / "task_executor")
    data = _source(n_features=40, n_hvg=35)
    data.uns.update(is_simulated=True, data_origin="synthetic_integrity_test")
    uri = "adata://native/normalized/v2"
    registry.register(uri, data, ArtifactType.ANNDATA, "native", "norm", "norm")
    integration_contract = TaskContract(task_id="integration", capability="integration", method="no_correction_v1",
                                        input_artifacts=[uri], parameters={"random_seed": 42, "n_components": 30})
    integrated = IntegrationCapability("no_correction_v1").execute(integration_contract, registry)
    capabilities = create_default_capability_registry()
    contract = TaskContract(task_id="clustering", capability="clustering", method="scanpy_leiden_umap_v1",
                            input_artifacts=integrated.output_artifacts, expected_outputs=["adata://native/annotated/v4"],
                            parameters={"random_seed": 13, "n_components": 20, "n_neighbors": 10})
    capabilities.prepare_contract(contract, registry)
    outcome = TaskExecutor(registry, capabilities).execute(contract)
    assert outcome.result.status.value == "success", outcome.result
    assert outcome.result.metrics["execution_attempts"] == 1
    signature = "native_clustering_integrity"
    registry.commit_task(outcome.staged, signature, outcome.result)
    outcome.staged.close()
    auditor = ScientificAuditor(additional_validators=capabilities.audit_validators())
    registry.begin_audit(signature, contract, outcome.result, auditor=auditor)
    report = auditor.audit_task(contract, outcome.result, registry)
    assert report.overall_passed, report
    registry.record_audit(signature, contract, outcome.result, report, auditor=auditor)
    reopened = ArtifactRegistry(registry.storage.base_dir)
    admitted = reopened.get_audited_many(signature=signature, contract=contract)
    assert len(admitted) == 1
    assert validate_clustering_execution_receipt(contract, outcome.result, reopened)["leiden"]["n_iterations"] == 2
    assert reopened.get_metadata(outcome.result.output_artifacts[0]).summary_metrics["is_simulated"] is True
    from eacbp.evidence.graph import EvidenceGraph
    from eacbp.report.markdown_report import ScientificReportGenerator
    from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest
    manifest = StudyManifest(study_id="native", biological_design=BiologicalDesign(species="human", tissue="test"),
                             data=DataSpec(raw_artifact_uri=uri))
    generator = ScientificReportGenerator(manifest, EvidenceGraph(), reopened, [outcome.result], [report])
    markdown = generator.generate_markdown()
    assert "### Actual clustering execution" in markdown
    assert "uncentered_sparse_truncated_svd" in markdown
    assert "| Actual representation dimensions | 30 |" in markdown
    assert "without refitting the graph" in markdown
    output_meta = reopened.get_metadata(outcome.result.output_artifacts[0])
    # File bytes and durable hashes are unchanged; receipt presentation must
    # still reject altered metadata parameters rather than trust the old pass.
    tampered_metadata = output_meta.model_copy(deep=True)
    tampered_metadata.parameters["resolution"] = 2.0
    reopened.registry[output_meta.uri] = tampered_metadata
    reopened._persist_index()
    with pytest.raises(ValueError, match="metadata parameters"):
        generator.generate_markdown()
    reopened.registry[output_meta.uri] = output_meta
    reopened._persist_index()
    import h5py
    with h5py.File(output_meta.storage_path, "r+") as artifact:
        artifact["uns/clustering/execution_receipt/leiden/resolution"][()] = 2.0
    with pytest.raises(Exception, match="integrity|hash|checksum"):
        generator.generate_markdown()


@pytest.mark.parametrize("attempts", [True, 0, 999, 1.5])
def test_untrusted_execution_attempts_are_rejected(native, attempts):
    registry, contract, result = native
    result = result.model_copy(deep=True)
    result.metrics["execution_attempts"] = attempts
    with pytest.raises(ValueError, match="Execution attempts"):
        validate_clustering_execution_receipt(contract, result, registry)


@pytest.mark.parametrize("corruption", ["metadata_resolution", "metadata_parents", "metadata_summary", "result_receipt", "source_sha"])
def test_cheap_receipt_reader_rejects_metadata_or_context_tampering(native, tmp_path, corruption):
    registry, contract, result = _rehashed(native, tmp_path, corruption)
    with pytest.raises((ValueError, AssertionError)):
        validate_clustering_execution_receipt(contract, result, registry, recompute_partition=False)
