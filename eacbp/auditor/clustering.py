"""Independently reconstruct the real Leiden execution from its input.

The third-party Scanpy implementation is reused; executor helpers are not.
This checks representation, graph, partition and provenance, not biological
truth or independent replication. UMAP is a display output outside this audit.
"""

from copy import deepcopy
import hashlib

import numpy as np
import pandas as pd
from scipy import sparse

from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity
from eacbp.capabilities.parameters import LeidenParameters
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.audit import model_to_jsonable, stable_fingerprint
from eacbp.schemas.task import TaskContract


def _digest(matrix):
    digest = hashlib.sha256()
    if sparse.issparse(matrix):
        value = sparse.csr_matrix(matrix, copy=True)
        value.sum_duplicates()
        value.sort_indices()
        value.eliminate_zeros()
        digest.update(stable_fingerprint({"shape": value.shape, "dtype": value.dtype.str, "format": "csr"}).encode())
        digest.update(np.ascontiguousarray(value.indptr.astype(np.int64)).tobytes())
        digest.update(np.ascontiguousarray(value.indices.astype(np.int64)).tobytes())
        digest.update(np.ascontiguousarray(value.data).tobytes())
    else:
        value = np.ascontiguousarray(matrix)
        digest.update(stable_fingerprint({"shape": value.shape, "dtype": value.dtype.str}).encode())
        digest.update(value.tobytes())
    return "sha256:" + digest.hexdigest()


def _partition(cell_ids, labels):
    # First occurrence determines a canonical code, independent of label names.
    mapping, membership = {}, []
    for label in labels:
        name = str(label)
        if name not in mapping:
            mapping[name] = len(mapping)
        membership.append(mapping[name])
    return stable_fingerprint({"cell_ids": list(map(str, cell_ids)), "partition": membership})


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _preserved_frame(before, after, columns=None):
    _require(before.index.astype(str).tolist() == after.index.astype(str).tolist(), "Cell or feature order changed")
    columns = list(before.columns) if columns is None else list(columns)
    _require(set(columns) <= set(after.columns), "Input metadata columns were removed")
    pd.testing.assert_frame_equal(before.loc[:, columns], after.loc[:, columns],
                                  check_dtype=False, check_categorical=False)


def _receipt_values(value):
    if isinstance(value, np.ndarray):
        return _receipt_values(value.tolist())
    if isinstance(value, dict):
        return {key: _receipt_values(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_receipt_values(item) for item in value]
    return model_to_jsonable(value)


def _inherited_receipt(data, metadata, registry):
    inherited = _receipt_values(data.uns.get("integration", {}).get("representation_receipt", {}))
    if not inherited:
        return {}
    _require(data.uns.get("integration", {}).get("method") == "no_correction_v1",
             "Unknown inherited representation-receipt method")
    _require(isinstance(inherited, dict), "Inherited representation receipt must be an object")
    source_uri = inherited.get("source_artifact_uri")
    _require(metadata.parent_uris == [source_uri], "Inherited representation lineage differs from its parent")
    source_meta, original = registry.get_scdata(source_uri)
    _require(_digest(original.X) == _digest(data.X), "No-correction input expression changed")
    _preserved_frame(original.var, data.var)
    _preserved_frame(original.obs, data.obs)
    mask = np.asarray(data.var["highly_variable"], dtype=bool) if "highly_variable" in data.var else np.ones(data.n_vars, dtype=bool)
    if "highly_variable" in data.var:
        _require(not data.var["highly_variable"].isna().any()
                 and data.var["highly_variable"].isin([True, False]).all(), "Invalid upstream feature mask")
    requested = metadata.parameters.get("n_components", 20)
    seed = metadata.parameters.get("random_seed", 42)
    is_sparse = sparse.issparse(original.X)
    actual = min(requested, max(1, data.n_obs - 2), max(1, int(mask.sum()) - 2)) if is_sparse else min(requested, data.n_obs, int(mask.sum()))
    algorithm = inherited.get("algorithm")
    routes = {
        "uncentered_sparse_truncated_svd": (False, "scipy.sparse.linalg.svds", False),
        "centered_dense_svd_fallback": (True, "numpy.linalg.svd", True),
        "centered_dense_svd": (True, "numpy.linalg.svd", False),
    }
    _require(algorithm in routes and ((is_sparse and algorithm != "centered_dense_svd")
                                     or (not is_sparse and algorithm == "centered_dense_svd")),
             "Inherited algorithm is inconsistent with its source matrix")
    centered, solver, fallback = routes[algorithm]
    expected = {
        "schema_version": 1, "method": "no_correction_v1", "source_artifact_uri": source_uri,
        "source_sha256": source_meta.sha256_hash,
        "cell_order_sha256": stable_fingerprint(original.obs.index.astype(str).tolist()),
        "feature_order_sha256": stable_fingerprint(original.var.index.astype(str).tolist()),
        "feature_selection": "highly_variable" if "highly_variable" in data.var else "all_features",
        "feature_mask_sha256": stable_fingerprint(mask.tolist()), "n_features_used": int(mask.sum()),
        "requested_n_components": int(requested), "actual_n_components": int(actual), "random_seed": int(seed),
        "representation_key": "X_pca", "representation_shape": list(data.obsm["X_pca"].shape),
        "representation_sha256": _digest(data.obsm["X_pca"]), "batch_correction_applied": False,
        "algorithm": algorithm, "zero_center": centered, "solver": solver, "fallback_used": fallback,
    }
    _require(data.obsm["X_pca"].shape == (data.n_obs, actual)
             and stable_fingerprint(inherited) == stable_fingerprint(expected),
             "Inherited representation settings do not match source axes, feature mask or embedding")
    return inherited


class ClusteringIntegrityValidator(BaseAuditor):
    def __init__(self):
        super().__init__("clustering_integrity_validator")

    def audit(self, contract, result, registry):
        report = ValidationReport(auditor_name=self.auditor_name, target_task_id=contract.task_id,
                                  target_artifact_uri=result.output_artifacts[0] if result.output_artifacts else None)
        if contract.capability != "clustering" or result.method_used != "scanpy_leiden_umap_v1":
            return report
        try:
            validate_clustering_execution_receipt(contract, result, registry, recompute_partition=True)
            report.add_check("leiden_execution_integrity", True, ValidationSeverity.ERROR,
                             "Reconstructed neighbors and Leiden partition from the input representation; checked effective settings, counts, axes and lineage. This is execution validation, not biological truth.",
                             {"neighbors_recomputed": True, "partition_recomputed": True, "umap_recomputed": False,
                              "upstream_representation_algorithm_recomputed": False,
                              "validation_scope": "execution_reconstruction_only", "independent_biological_validation": False})
        except Exception as exc:
            report.add_check("leiden_execution_integrity", False, ValidationSeverity.ERROR,
                             f"Independent Leiden reconstruction failed: {exc}",
                             {"validation_scope": "execution_reconstruction_only", "independent_biological_validation": False})
        return report

    def _verify(self, contract, result, registry, *, recompute_partition):
        if recompute_partition:
            import scanpy as sc

        _require(contract.method in (None, "scanpy_leiden_umap_v1")
                 and len(contract.input_artifacts) == len(result.output_artifacts) == 1,
                 "Real Leiden requires one input and one output under its actual method")
        _require(result.input_artifacts == contract.input_artifacts and result.task_id == contract.task_id
                 and result.capability == "clustering" and result.status.value == "success"
                 and (not contract.expected_outputs or result.output_artifacts == contract.expected_outputs)
                 and result.inference_contract is None and result.scientific_result is None,
                 "Result input or task identity differs from contract")
        source_meta, source = registry.get_scdata(contract.input_artifacts[0])
        output_meta, output = registry.get_scdata(result.output_artifacts[0])
        _require(output_meta.type == ArtifactType.ANNDATA
                 and output_meta.created_by_task == contract.task_id
                 and output_meta.operation == "scanpy_leiden_umap_v1"
                 and output_meta.parent_uris == contract.input_artifacts
                 and stable_fingerprint(output_meta.parameters) == stable_fingerprint(contract.parameters),
                 "Output metadata parameters, method or lineage differ from contract")
        params = contract.parameters
        LeidenParameters.model_validate(params)
        resolution = params.get("resolution", 1)
        _require(isinstance(resolution, (int, float)) and not isinstance(resolution, bool)
                 and np.isfinite(resolution) and resolution > 0, "resolution must be finite and positive")
        seed = int(params.get("random_seed", 42))
        _require(output_meta.random_seed == seed, "Metadata random seed differs from actual Leiden seed")
        _require(min(source.shape) >= 3 and source.shape == output.shape
                 and source.obs.index.is_unique and source.var.index.is_unique,
                 "Invalid source dimensions or ambiguous cell/feature identities")
        _require(source.obs.index.astype(str).is_unique and source.var.index.astype(str).is_unique,
                 "Source identities collide after conversion to strings")
        _require(_digest(source.X) == _digest(output.X), "Input expression values changed during clustering")
        _preserved_frame(source.var, output.var)
        _preserved_frame(source.obs, output.obs, [c for c in source.obs if c not in {"cluster", "leiden", "cell_type"}])
        _require(set(source.layers) == set(output.layers), "Clustering changed the input layer set")
        for layer in source.layers:
            _require(_digest(source.layers[layer]) == _digest(output.layers[layer]), f"Input layer {layer} changed")
        for name in source.obsm:
            if name not in {"X_umap", "X_embedding_2d"}:
                _require(name in output.obsm and _digest(source.obsm[name]) == _digest(output.obsm[name]),
                         f"Input embedding {name} changed")
        for name in source.obsp:
            if name not in {"connectivities", "distances"}:
                _require(name in output.obsp and _digest(source.obsp[name]) == _digest(output.obsp[name]),
                         f"Input pairwise matrix {name} changed")

        rebuilt = source.to_anndata()
        reused = "X_pca" in source.obsm
        if not reused:
            components = min(params.get("n_components", 20), source.n_obs - 1, source.n_vars - 1)
            if recompute_partition:
                sc.pp.pca(rebuilt, n_comps=components, zero_center=True, svd_solver="arpack", random_state=seed, mask_var=None)
            else:
                _require(output.obsm["X_pca"].shape == (source.n_obs, components), "Computed PCA shape differs from contract")
                rebuilt.obsm["X_pca"] = output.obsm["X_pca"].copy()
        _require(np.isfinite(rebuilt.obsm["X_pca"]).all()
                 and _digest(rebuilt.obsm["X_pca"]) == _digest(output.obsm["X_pca"]),
                 "Actual clustering representation differs from the input or independent PCA")
        neighbor_count = min(params.get("n_neighbors", 15), source.n_obs - 1)
        expected_neighbor_settings = {
            "connectivities_key": "connectivities", "distances_key": "distances",
            "params": {"n_neighbors": neighbor_count, "method": "umap", "random_state": seed,
                       "metric": "euclidean", "use_rep": "X_pca"},
        }
        if recompute_partition:
            sc.pp.neighbors(rebuilt, n_neighbors=neighbor_count, use_rep="X_pca", random_state=seed)
        else:
            rebuilt.uns["neighbors"] = expected_neighbor_settings
            for name in ("connectivities", "distances"):
                rebuilt.obsp[name] = output.obsp[name].copy()
        for name in ("connectivities", "distances"):
            _require(name in output.obsp and _digest(rebuilt.obsp[name]) == _digest(output.obsp[name]),
                     f"Neighbor {name} differs from independent reconstruction")
        _require(stable_fingerprint(rebuilt.uns["neighbors"]) == stable_fingerprint(output.uns.get("neighbors", {})),
                 "Persisted Scanpy neighbor settings differ from actual settings")
        _require(stable_fingerprint(rebuilt.uns["neighbors"]) == stable_fingerprint(expected_neighbor_settings),
                 "Actual Scanpy neighbor settings differ from the fixed implementation contract")
        if recompute_partition:
            sc.tl.leiden(rebuilt, resolution=float(resolution), random_state=seed, flavor="leidenalg",
                         key_added="cluster", directed=False, n_iterations=2)
        _require("cluster" in output.obs and "leiden" in output.obs
                 and not output.obs["cluster"].isna().any() and not output.obs["leiden"].isna().any()
                 and output.obs["cluster"].nunique() == output.obs["cluster"].astype(str).nunique()
                 and output.obs["cluster"].astype(str).tolist() == output.obs["leiden"].astype(str).tolist(),
                 "Cluster labels are missing, incomplete or disagree with the Leiden alias")
        partition = _partition(output.obs.index, output.obs["cluster"])
        if recompute_partition:
            _require(partition == _partition(rebuilt.obs_names, rebuilt.obs["cluster"]),
                     "Cell partition differs from independent Leiden reconstruction")
        else:
            rebuilt.obs["cluster"] = output.obs["cluster"].copy()
        _require(output.obs["cluster"].nunique() == rebuilt.obs["cluster"].nunique(),
                 "Reported cluster count differs from the actual partition")
        _require(stable_fingerprint(output.uns.get("cluster", {}).get("params", {}))
                 == stable_fingerprint({"resolution": float(resolution), "random_state": seed, "n_iterations": 2}),
                 "Persisted Scanpy Leiden settings differ from actual settings")

        markers = params.get("marker_reference", {})
        _require(isinstance(markers, dict) and all(isinstance(key, str) and key.strip()
                 and isinstance(values, list) and values and all(isinstance(g, str) and g.strip() for g in values)
                 for key, values in markers.items()), "Invalid supplied marker reference")
        if markers:
            genes = output.var["gene_name"].astype(str).to_numpy() if "gene_name" in output.var else output.var.index.astype(str).to_numpy()
            expected_labels = {}
            for cluster in output.obs["cluster"].unique():
                inside = (output.obs["cluster"] == cluster).to_numpy()
                means = np.asarray(source.X[inside].mean(axis=0)).ravel()
                scores = {name: float(means[np.isin(genes, reference)].mean()) for name, reference in markers.items()
                          if np.isin(genes, reference).any()}
                expected_labels[cluster] = max(scores, key=scores.get) if scores and max(scores.values()) > 0 else f"Cluster_{cluster}"
            expected_types = [expected_labels[c] for c in output.obs["cluster"]]
            _require(output.obs["cell_type"].astype(str).tolist() == expected_types, "Marker annotations were not reconstructed")
        elif "cell_type" in source.obs:
            _preserved_frame(source.obs, output.obs, ["cell_type"])
        else:
            _require(output.obs["cell_type"].astype(str).tolist() == [f"Cluster_{c}" for c in output.obs["cluster"]],
                     "Unannotated cluster labels were replaced by invented identities")

        inherited = _inherited_receipt(source, source_meta, registry) if reused else {}
        integration_method = str(source.uns.get("integration", {}).get("method", "unavailable"))
        role_columns = [name for name in ("donor", "donor_id", "batch", "batch_id") if name in source.obs]
        role_values = {name: [None if pd.isna(value) else str(value) for value in source.obs[name]] for name in role_columns}
        annotation = "user_marker_reference" if markers else "existing_obs_or_unannotated_clusters"
        expected = {
            "schema_version": 1, "schema_name": "leiden_execution_receipt_v1", "method": "scanpy_leiden_umap_v1",
            "source_artifact_uri": contract.input_artifacts[0], "source_sha256": source_meta.sha256_hash,
            "cell_order_sha256": stable_fingerprint(source.obs.index.astype(str).tolist()),
            "feature_order_sha256": stable_fingerprint(source.var.index.astype(str).tolist()),
            "counts_source": "layers.counts" if "counts" in source.layers else "X_input_representation",
            "counts_sha256": _digest(source.layers["counts"] if "counts" in source.layers else source.X),
            "representation": {
                "key": "X_pca", "source": "input.obsm.X_pca" if reused else "computed_from_input.X",
                "reused": reused, "source_integration_method": integration_method,
                "requested_n_components": int(params.get("n_components", 20)),
                "actual_n_components": int(rebuilt.obsm["X_pca"].shape[1]),
                "shape": list(rebuilt.obsm["X_pca"].shape), "sha256": _digest(rebuilt.obsm["X_pca"]),
                "algorithm": inherited.get("algorithm", "inherited_not_documented") if reused else "scanpy_centered_pca_arpack",
                "feature_selection": inherited.get("feature_selection", "inherited_not_documented") if reused else "all_features",
                "feature_mask_sha256": inherited.get("feature_mask_sha256", "") if reused else stable_fingerprint([True] * source.n_vars),
                "inherited_representation_receipt": inherited,
            },
            "neighbors": {
                "requested_n_neighbors": int(params.get("n_neighbors", 15)),
                "actual_n_neighbors": int(rebuilt.uns["neighbors"]["params"]["n_neighbors"]),
                "use_rep": "X_pca", "actual_n_components": int(rebuilt.obsm["X_pca"].shape[1]),
                "all_representation_components_used": True, "metric": "euclidean", "method": "umap", "random_state": seed,
                "connectivities_sha256": _digest(rebuilt.obsp["connectivities"]), "distances_sha256": _digest(rebuilt.obsp["distances"]),
            },
            "leiden": {
                "resolution": float(resolution), "random_state": seed, "flavor": "leidenalg",
                "partition_type": "RBConfigurationVertexPartition", "directed": False, "use_weights": True, "n_iterations": 2,
                "n_clusters": int(rebuilt.obs["cluster"].nunique()), "partition_sha256": partition,
            },
            "input_roles": {"columns": role_columns, "metadata_sha256": stable_fingerprint(role_values), "donor_independence_established": False},
            "annotation_source": annotation, "clustering_reference_labels_used": False,
            "independent_biological_validation": False, "validation_scope": "execution_reconstruction_only", "umap_reconstruction_in_scope": False,
        }
        summary = {"n_cells": source.n_obs, "n_clusters": int(rebuilt.obs["cluster"].nunique()),
                   "annotation_source": annotation, "execution_receipt": expected}
        actual = _receipt_values(output.uns.get("clustering", {}))
        _require(stable_fingerprint(actual) == stable_fingerprint({"method": "scanpy_leiden_umap_v1", **summary}),
                 "Payload execution receipt or clustering summary differs from independent reconstruction")
        inherited_provenance = {key: deepcopy(source_meta.summary_metrics[key])
                                for key in ("is_simulated", "data_origin") if key in source_meta.summary_metrics}
        _require(stable_fingerprint(output_meta.summary_metrics) == stable_fingerprint({**summary, **inherited_provenance}),
                 "Metadata summary metrics differ from the verified payload receipt")
        routing = {key: params[key] for key in ("target_cell_type", "target_branch", "fdr_family") if key in params}
        result_metrics = dict(result.metrics)
        if "execution_attempts" in result_metrics:
            attempts = result_metrics.pop("execution_attempts")
            _require(isinstance(attempts, int) and not isinstance(attempts, bool)
                     and 1 <= attempts <= contract.retry_policy.require_human_after,
                     "Execution attempts are outside the contracted retry limit")
        _require(stable_fingerprint(result_metrics) in {stable_fingerprint(summary), stable_fingerprint({**summary, **routing})},
                 "Result metrics differ from the verified payload receipt")
        return deepcopy(expected)


def validate_clustering_execution_receipt(contract, result, registry, *, recompute_partition=False):
    """Validate settings and lineage, optionally reconstructing graph/partition.

The cheap path is for presenting already durably admitted outputs. It binds
their persisted graph and partition hashes without claiming a new graph fit.
The admission auditor always uses ``recompute_partition=True``.
"""
    if isinstance(contract, dict):
        contract = TaskContract.model_validate(contract)
    _require(contract.capability == "clustering" and result.method_used == "scanpy_leiden_umap_v1",
             "A real Leiden contract/result is required")
    return ClusteringIntegrityValidator()._verify(contract, result, registry, recompute_partition=recompute_partition)
