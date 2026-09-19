"""Real, optional-library implementations. No substitution with baseline algorithms."""
from importlib.metadata import version
import numpy as np
import pandas as pd
from scipy import sparse, stats
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.numerics import benjamini_hochberg
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskResult, TaskStatus
from eacbp.artifact.uri import ArtifactURI


def _input(contract, registry):
    _, payload = registry.get(contract.input_artifacts[0])
    data = SCData.from_dict(payload)
    if min(data.shape) < 3:
        raise ValueError("Standard methods require at least three cells and genes")
    values = data.X.data if sparse.issparse(data.X) else data.X
    if not np.isfinite(values).all():
        raise ValueError("Expression must be finite")
    return data.to_anndata()


def _pca(adata, params):
    import scanpy as sc
    if "X_pca" not in adata.obsm:
        n = min(int(params.get("n_components", 20)), adata.n_obs - 1, adata.n_vars - 1)
        sc.pp.pca(adata, n_comps=n, zero_center=True, svd_solver="arpack",
                  random_state=int(params.get("random_seed", 42)), mask_var=None)
    if not np.isfinite(adata.obsm["X_pca"]).all():
        raise ValueError("PCA produced a nonfinite embedding")


def _neighbors(adata, params):
    import scanpy as sc
    _pca(adata, params)
    sc.pp.neighbors(adata, n_neighbors=min(int(params.get("n_neighbors", 15)), adata.n_obs - 1),
                    use_rep="X_pca", random_state=int(params.get("random_seed", 42)))


class _LibraryCapability(BaseCapability):
    def _save(self, contract, registry, payload, default_path, metrics, dependencies, artifact_type=ArtifactType.ANNDATA):
        sid = ArtifactURI.parse(contract.input_artifacts[0]).study_id
        uri = contract.expected_outputs[0] if contract.expected_outputs else f"{'table' if artifact_type == ArtifactType.TABLE else 'adata'}://{sid}/{default_path}"
        registry.register(uri, payload, artifact_type, sid, contract.task_id, self.implementation_id,
                          parent_uris=contract.input_artifacts, parameters=contract.parameters,
                          software_versions={name: version(name) for name in dependencies},
                          random_seed=int(contract.parameters.get("random_seed", 42)), summary_metrics=metrics)
        return TaskResult(task_id=contract.task_id, capability=self.capability_name, status=TaskStatus.SUCCESS,
                          method_used=self.implementation_id, input_artifacts=contract.input_artifacts,
                          output_artifacts=[uri], executed_operations=self.contract_operations, metrics=metrics)


class HarmonyIntegrationCapability(_LibraryCapability):
    contract_operations = ["compute_pca", "harmony_batch_integration", "evaluate_batch_mixing"]

    def __init__(self):
        super().__init__("integration", "harmonypy_v1", ImplementationType.PYTHON_TOOL,
                         accepts_types=[ArtifactType.ANNDATA], output_types=[ArtifactType.ANNDATA])
        self.legacy_aliases = {"harmony": "harmonypy_v1", "harmony_v1": "harmonypy_v1"}

    def execute(self, contract, registry):
        import harmonypy
        from eacbp.capabilities.integration import calculate_batch_mixing_score
        adata = _input(contract, registry)
        params = contract.parameters
        batch_col = params.get("batch_col", "batch")
        if batch_col not in adata.obs or adata.obs[batch_col].isna().any():
            raise ValueError("Harmony requires complete batch metadata")
        if adata.obs[batch_col].nunique() < 2:
            raise ValueError("Harmony requires at least two batches; select no_correction_v1")
        _pca(adata, params)
        before = np.asarray(adata.obsm["X_pca"]).copy()
        seed = int(params.get("random_seed", 42))
        harmony = harmonypy.run_harmony(before, adata.obs, batch_col,
                    nclust=min(max(2, adata.n_obs // 30), adata.n_obs - 1),
                    random_state=seed, max_iter_harmony=int(params.get("max_iter_harmony", 10)), verbose=False)
        corrected = np.asarray(harmony.Z_corr)
        # harmonypy 0.0.10 returns components x cells; this version is pinned.
        corrected = corrected.T
        if corrected.shape != before.shape or not np.isfinite(corrected).all():
            raise ValueError("Harmony returned invalid embedding")
        adata.obsm["X_pca_uncorrected"] = before
        adata.obsm["X_pca"] = corrected
        adata.obsm["X_pca_harmony"] = corrected.copy()
        batches = adata.obs[batch_col].astype(str).to_numpy()
        metrics = {"pre_batch_mixing": calculate_batch_mixing_score(before, batches, random_seed=seed),
                   "post_batch_mixing": calculate_batch_mixing_score(corrected, batches, random_seed=seed),
                   "batch_correction_applied": True, "n_cells": adata.n_obs}
        adata.uns["integration"] = {"method": self.implementation_id, "batch_col": batch_col, **metrics}
        return self._save(contract, registry, SCData.from_anndata(adata), "integrated/v3", metrics, ["harmonypy", "scanpy", "anndata"])


class LeidenClusteringCapability(_LibraryCapability):
    contract_operations = ["compute_pca", "build_neighbor_graph", "leiden_clustering", "umap_embedding", "marker_score_annotation"]

    def __init__(self):
        super().__init__("clustering", "scanpy_leiden_umap_v1", ImplementationType.PYTHON_TOOL,
                         accepts_types=[ArtifactType.ANNDATA], output_types=[ArtifactType.ANNDATA])
        self.legacy_aliases = {"leiden_knn_v1": "scanpy_leiden_umap_v1"}

    def execute(self, contract, registry):
        import scanpy as sc
        adata = _input(contract, registry)
        params = contract.parameters
        seed = int(params.get("random_seed", 42))
        _neighbors(adata, params)
        sc.tl.leiden(adata, resolution=float(params.get("resolution", 1)), random_state=seed,
                     flavor="leidenalg", key_added="cluster", directed=False, n_iterations=2)
        sc.tl.umap(adata, random_state=seed)
        adata.obs["leiden"] = adata.obs["cluster"].copy()
        adata.obsm["X_embedding_2d"] = adata.obsm["X_umap"].copy()
        # Explicit marker reference only; existing measured annotations may be retained.
        markers = params.get("marker_reference", {})
        if markers:
            genes = adata.var["gene_name"].astype(str).to_numpy() if "gene_name" in adata.var else adata.var_names.to_numpy()
            labels = {}
            for cluster in adata.obs["cluster"].unique():
                mask = (adata.obs["cluster"] == cluster).to_numpy()
                means = np.asarray(adata.X[mask].mean(axis=0)).ravel()
                scores = {name: float(means[np.isin(genes, reference)].mean())
                          for name, reference in markers.items() if np.isin(genes, reference).any()}
                labels[cluster] = max(scores, key=scores.get) if scores and max(scores.values()) > 0 else f"Cluster_{cluster}"
            adata.obs["cell_type"] = [labels[c] for c in adata.obs["cluster"]]
        elif "cell_type" not in adata.obs:
            adata.obs["cell_type"] = [f"Cluster_{c}" for c in adata.obs["cluster"]]
        metrics = {"n_cells": adata.n_obs, "n_clusters": int(adata.obs["cluster"].nunique()),
                   "annotation_source": "user_marker_reference" if markers else "existing_obs_or_unannotated_clusters"}
        adata.uns["clustering"] = {"method": self.implementation_id, **metrics}
        return self._save(contract, registry, SCData.from_anndata(adata), "annotated/v4", metrics, ["scanpy", "leidenalg", "igraph", "umap-learn", "anndata"])


class DPTTrajectoryCapability(_LibraryCapability):
    contract_operations = ["compute_pca", "build_neighbor_graph", "diffusion_map", "diffusion_pseudotime", "paga_graph", "spearman_gene_association", "benjamini_hochberg_correction"]

    def __init__(self):
        super().__init__("trajectory_inference", "scanpy_dpt_v1", ImplementationType.PYTHON_TOOL,
                         accepts_types=[ArtifactType.ANNDATA], output_types=[ArtifactType.TABLE])

    def execute(self, contract, registry):
        import scanpy as sc
        adata = _input(contract, registry)
        params = contract.parameters
        root = params.get("root_cell_id")
        ids = adata.obs["cell_id"].astype(str).to_numpy() if "cell_id" in adata.obs else adata.obs_names.to_numpy()
        hits = np.flatnonzero(ids == str(root)) if root is not None else []
        if len(hits) != 1:
            raise ValueError("DPT requires root_cell_id identifying exactly one observed cell")
        _neighbors(adata, params)
        if sparse.csgraph.connected_components(adata.obsp["connectivities"], directed=False)[0] != 1:
            raise ValueError("DPT graph is disconnected; analyze a connected population")
        n_dcs = min(int(params.get("n_dcs", 10)), adata.n_obs - 1)
        sc.tl.diffmap(adata, n_comps=n_dcs, random_state=int(params.get("random_seed", 42)))
        adata.uns["iroot"] = int(hits[0])
        sc.tl.dpt(adata, n_dcs=n_dcs)
        if params.get("run_paga", False):
            group_key = params.get("paga_groups", "cluster")
            if group_key not in adata.obs or adata.obs[group_key].nunique() < 2:
                raise ValueError("PAGA requires an observed grouping with at least two groups")
            adata.obs[group_key] = adata.obs[group_key].astype("category")
            sc.tl.paga(adata, groups=group_key)
        pseudotime = adata.obs["dpt_pseudotime"].to_numpy()
        if not np.isfinite(pseudotime).all():
            raise ValueError("DPT returned nonfinite pseudotime")
        rows = []
        genes = adata.var["gene_name"].astype(str).to_numpy() if "gene_name" in adata.var else adata.var_names.to_numpy()
        for j, gene in enumerate(genes):
            values = adata.X[:, j]
            values = values.toarray().ravel() if sparse.issparse(values) else np.asarray(values).ravel()
            rho, p = stats.spearmanr(pseudotime, values) if np.ptp(values) > 0 else (0., 1.)
            rows.append({"gene": gene, "spearman_rho": float(rho) if np.isfinite(rho) else 0., "p_value": float(p) if np.isfinite(p) else 1.})
        table = pd.DataFrame(rows)
        table["fdr_q_value"] = benjamini_hochberg(table["p_value"].to_numpy())
        metrics = {"root_cell_id": str(root), "n_cells": adata.n_obs, "root_selection": "user_supplied",
                   "stability_evaluated": False}
        result = self._save(contract, registry, table, "trajectory_results/v1", metrics, ["scanpy", "anndata"], ArtifactType.TABLE)
        if not params.get("run_paga", False):
            result.executed_operations = [op for op in result.executed_operations if op != "paga_graph"]
        sid = ArtifactURI.parse(contract.input_artifacts[0]).study_id
        # CellRank plans append the persisted pseudotime AnnData output to the
        # trajectory contract.  Honor that namespace so target branches do
        # not overwrite one another; retain the legacy URI as fallback.
        uri = (
            contract.expected_outputs[1]
            if len(contract.expected_outputs) > 1
            else f"adata://{sid}/diffusion_pseudotime/v1"
        )
        registry.register(uri, SCData.from_anndata(adata), ArtifactType.ANNDATA, sid, contract.task_id,
                          self.implementation_id, parent_uris=contract.input_artifacts, parameters=params,
                          summary_metrics=metrics, software_versions={"scanpy": version("scanpy")})
        result.output_artifacts.append(uri)
        return result


class CellRankFateCapability(_LibraryCapability):
    """CellRank absorption probabilities with explicitly supplied terminal cells."""
    contract_operations = ["cellrank_transition_kernel", "cellrank_set_terminal_states", "cellrank_fate_probabilities"]

    def __init__(self):
        super().__init__("fate_mapping", "cellrank_fate_v1", ImplementationType.PYTHON_TOOL,
                         accepts_types=[ArtifactType.ANNDATA], output_types=[ArtifactType.TABLE])

    def execute(self, contract, registry):
        import cellrank as cr
        adata = _input(contract, registry)
        params = contract.parameters
        terminals = params.get("terminal_states")
        if not isinstance(terminals, dict) or not terminals or any(not ids for ids in terminals.values()):
            raise ValueError("CellRank requires terminal_states mapping fate names to observed cell IDs")
        if "cell_id" in adata.obs:
            adata.obs_names = adata.obs["cell_id"].astype(str)
        if not adata.obs_names.is_unique:
            raise ValueError("CellRank requires unique cell IDs")
        selected = [str(cell) for group in terminals.values() for cell in group]
        if len(selected) != len(set(selected)) or not set(selected) <= set(adata.obs_names):
            raise ValueError("Terminal cells must be distinct observed cells")
        kernel_kind = params.get("kernel", "pseudotime")
        if kernel_kind == "precomputed":
            key = params.get("transition_key", "transition_matrix")
            if key not in adata.obsp:
                raise ValueError("Precomputed CellRank requires an observed transition matrix in obsp")
            matrix = sparse.csr_matrix(adata.obsp[key])
            if not np.isfinite(matrix.data).all() or (matrix.data < 0).any() or not np.allclose(np.asarray(matrix.sum(axis=1)).ravel(), 1):
                raise ValueError("Transition matrix must be finite, nonnegative and row-stochastic")
            kernel = cr.kernels.PrecomputedKernel(matrix, adata=adata)
        elif kernel_kind == "pseudotime":
            key = params.get("time_key", "dpt_pseudotime")
            if key not in adata.obs or not np.isfinite(adata.obs[key].to_numpy()).all():
                raise ValueError("CellRank pseudotime kernel requires finite observed pseudotime")
            if "connectivities" not in adata.obsp:
                raise ValueError("CellRank pseudotime kernel requires a neighbor graph")
            kernel = cr.kernels.PseudotimeKernel(adata, time_key=key)
            kernel.compute_transition_matrix(threshold_scheme="soft", n_jobs=1, show_progress_bar=False)
        else:
            raise ValueError("kernel must be pseudotime or precomputed; raw velocity inference is not substituted")
        estimator = cr.estimators.GPCCA(kernel)
        estimator.set_terminal_states({str(name): list(map(str, ids)) for name, ids in terminals.items()})
        estimator.compute_fate_probabilities(solver="gmres", use_petsc=False, n_jobs=1, show_progress_bar=False)
        probabilities = np.asarray(estimator.fate_probabilities)
        if not np.isfinite(probabilities).all() or (probabilities < -1e-6).any() or not np.allclose(probabilities.sum(axis=1), 1, atol=1e-3):
            raise ValueError("CellRank returned invalid fate probabilities")
        table = pd.DataFrame(probabilities, index=adata.obs_names, columns=estimator.fate_probabilities.names)
        metrics = {"n_cells": adata.n_obs, "n_fates": table.shape[1], "kernel": kernel_kind,
                   "terminal_state_source": "user_supplied", "model_conditional_probabilities": True}
        return self._save(contract, registry, table, "fate_probabilities/v1", metrics, ["cellrank", "pygpcca", "anndata"], ArtifactType.TABLE)
