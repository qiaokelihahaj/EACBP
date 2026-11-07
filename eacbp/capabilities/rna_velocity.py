"""Optional native scVelo stochastic RNA-velocity analysis.

Spliced/unspliced counts, project-generated neighborhoods, and kinetic model
inference are the inputs. Published roots, cell types, and measured time labels
are never supplied to the velocity model. Model direction and pseudotime are
within-dataset descriptions, not measured time or causal fate observations.
"""

from __future__ import annotations

import importlib
from importlib import metadata
import re
from typing import Literal

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, model_validator
from scipy import sparse

from eacbp.artifact.uri import ArtifactURI
from eacbp.capabilities.base import BaseCapability, CapabilityDescriptor
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import EvidenceNode, EvidenceStrength, EvidenceType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


MODEL_ASSUMPTIONS = [
    "Spliced and unspliced layers represent aligned nonnegative integer RNA counts.",
    "A stochastic transcription/splicing/degradation model describes the observed RNA kinetics.",
    "Local project-generated transcriptomic neighborhoods support first and second RNA moments.",
    "Kinetic parameters and inferred direction can be biased by model mismatch and sampling.",
    "Model-inferred velocity and pseudotime are not measured time, observed fate, or causal evidence.",
]


class RNAVelocityParameters(BaseModel):
    model_config = ConfigDict(extra="allow")

    groupby: Literal["cluster"] = "cluster"
    n_jobs: Literal[1] = 1
    n_dcs: int = Field(10, ge=2, le=100)
    target_sum: float = Field(10000.0, gt=0, allow_inf_nan=False)
    random_seed: int = Field(42, ge=0)

    @model_validator(mode="before")
    @classmethod
    def reject_external_orientation(cls, values):
        forbidden = {"root", "root_key", "root_cell_id", "iroot", "end_key", "end_cell_id", "latent_time",
                     "mode", "filter_genes", "use_highly_variable", "spliced_layer", "unspliced_layer",
                     "fit_offset", "fit_offset2", "min_r2", "min_likelihood", "perc", "use_raw",
                     "groups", "tkey", "gene_subset", "vkey", "n_neighbors", "n_pcs", "basis", "approx"}
        if isinstance(values, dict) and forbidden.intersection(values):
            raise ValueError("RNA velocity fixes its stochastic model and infers orientation from kinetics; caller roots, layer substitutions and unimplemented scientific overrides are unsupported")
        return values


def _raw_counts(matrix, name):
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.rint(values)).all():
        raise ValueError(f"RNA velocity requires finite nonnegative integer {name} counts")
    if not np.any(values > 0):
        raise ValueError(f"RNA velocity requires nonempty {name} counts")


def _require_velocity_numpy():
    match = re.match(r"(\d+)\.(\d+)", np.__version__)
    if match is None or tuple(map(int, match.groups())) >= (2, 4):
        raise ImportError(
            f"scVelo 0.3.x stochastic fitting is incompatible with NumPy {np.__version__}; "
            "it requires numpy>=1.24,<2.4. Install EACBP[velocity] in a compatible environment. "
            "No substitute RNA-velocity method is used."
        )


def _velocity_plan(*, manifest, parameters, tasks):
    if any(task.capability == "rna_velocity" for task in tasks):
        return []
    clustering = next((task for task in reversed(tasks) if task.capability == "clustering"), None)
    if clustering is None or clustering.method not in {None, "scanpy_leiden_umap_v1"}:
        raise ValueError("RNA velocity requires native Scanpy clustering with project-generated neighbors/PCA")
    return TaskContract(
        task_id="task_005c_rna_velocity", capability="rna_velocity", method="scvelo_stochastic_v1",
        depends_on=[clustering.task_id], input_artifacts=[clustering.expected_outputs[0]], parameters=parameters,
        allowed_operations=list(RNAVelocityCapability.contract_operations),
        forbidden_operations=["filter_cells", "filter_genes", "recluster", "replace_cell_type_labels", "recover_dynamics", "latent_time"],
        expected_outputs=[f"adata://{manifest.study_id}/rna_velocity/v1", f"table://{manifest.study_id}/rna_velocity_summary/v1"],
        validation_requirements=["rna_velocity_integrity"],
    )


def _velocity_evidence(contract, result, report, registry):
    data = registry.get_scdata(result.output_artifacts[0])[1]
    return [EvidenceNode(
        evidence_id=f"E_{contract.task_id}_rna_velocity", type=EvidenceType.RNA_VELOCITY,
        strength=EvidenceStrength.WEAK, score=0.5, audit_passed=True,
        summary=("scVelo stochastic kinetics inferred local RNA velocity, model root/end scores and velocity pseudotime. "
                 "Direction depends on kinetic and neighborhood assumptions; it is not measured time, observed fate, or causal evidence."),
        source_task_id=contract.task_id, source_artifact_uris=list(result.output_artifacts),
        is_simulated=bool(data.uns.get("is_simulated", False)), validation_scope="within_dataset_model_inferred",
        data_origin_uris=list(contract.input_artifacts),
        metrics={**result.metrics, "independent_replication": False, "confirmatory_allowed": False},
    )]


class RNAVelocityCapability(BaseCapability):
    contract_operations = ["preserve_rna_counts", "normalize_rna_layers", "compute_rna_moments",
                           "stochastic_rna_velocity", "velocity_graph", "infer_kinetic_terminal_states",
                           "velocity_pseudotime", "velocity_confidence"]

    def __init__(self):
        from eacbp.auditor.rna_velocity import RNAVelocityValidator

        super().__init__(
            "rna_velocity", "scvelo_stochastic_v1", accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA, ArtifactType.TABLE],
            descriptor=CapabilityDescriptor(
                "rna_velocity", "scvelo_stochastic_v1", parameter_model=RNAVelocityParameters,
                input_types=[ArtifactType.ANNDATA], output_types=[ArtifactType.ANNDATA, ArtifactType.TABLE],
                required_audit_ids=["rna_velocity_integrity"], validator=RNAVelocityValidator(),
                evidence_extractor=_velocity_evidence, plan_factory=_velocity_plan, scope="shared",
                runtime_packages=("anndata", "scvelo"),
                description="Native scVelo stochastic velocity with preserved RNA counts and model-inferred orientation.",
            ),
        )

    def execute(self, contract, registry):
        parameters = RNAVelocityParameters.model_validate(contract.parameters)
        if len(contract.input_artifacts) != 1:
            raise ValueError("RNA velocity requires exactly one native annotated single-cell input")
        in_uri = contract.input_artifacts[0]
        data = registry.get_scdata(in_uri)[1]
        if data.n_obs < 4 or data.n_vars < 3:
            raise ValueError("RNA velocity requires at least four cells and three genes")
        if data.uns.get("clustering", {}).get("method") != "scanpy_leiden_umap_v1":
            raise ValueError("RNA velocity requires project-generated native Scanpy clustering")
        if "cluster" not in data.obs or data.obs["cluster"].isna().any():
            raise ValueError("RNA velocity requires complete project cluster labels")
        if "X_pca" not in data.obsm or not np.isfinite(data.obsm["X_pca"]).all() or "neighbors" not in data.uns:
            raise ValueError("RNA velocity requires the project's finite PCA and neighbor metadata")
        for name in ("connectivities", "distances"):
            graph = data.obsp.get(name)
            if graph is None or graph.shape != (data.n_obs, data.n_obs):
                raise ValueError(f"RNA velocity requires the project's aligned {name} graph")
            values = graph.data if sparse.issparse(graph) else np.asarray(graph)
            if not np.isfinite(values).all() or np.any(values < 0):
                raise ValueError(f"Invalid project neighbor {name} values")
        for name in ("counts", "spliced", "unspliced"):
            if name not in data.layers:
                raise ValueError(f"RNA velocity requires an aligned raw {name} layer")
            _raw_counts(data.layers[name], name)
        _require_velocity_numpy()
        try:
            scv = importlib.import_module("scvelo")
        except ImportError as exc:
            raise ImportError("Native RNA velocity requires the velocity extra (scvelo>=0.3.3,<0.4); no substitute method is used") from exc
        adata = data.to_anndata()
        original_x = adata.X.copy()
        # scVelo normalizes X as well as RNA layers and uses gene-count
        # correlation in root diagnostics; give that step actual raw counts.
        adata.X = data.layers["counts"].copy()
        for name in list(adata.obs):
            if name in {"root_cells", "end_points", "latent_time", "initial_size", "initial_size_spliced", "initial_size_unspliced",
                        "cell_fate", "cell_fate_confidence", "cell_origin", "cell_origin_confidence"} or name.startswith("velocity_"):
                del adata.obs[name]
        for name in list(adata.var):
            if name in {"gene_count_corr", "spearmans_score"} or name.startswith(("velocity_", "fit_")):
                del adata.var[name]
        for name in list(adata.uns):
            if name == "iroot" or name.startswith(("velocity_", "recover_dynamics", "rank_velocity")):
                del adata.uns[name]
        for name in ("Ms", "Mu", "Mss", "Mus", "velocity", "velocity_u", "variance_velocity"):
            if name in adata.layers:
                del adata.layers[name]
        adata.layers["spliced_counts"] = data.layers["spliced"].copy()
        adata.layers["unspliced_counts"] = data.layers["unspliced"].copy()
        adata.obs["cluster"] = adata.obs["cluster"].astype("category")
        # All stochastic routines share a reproducible seed, restored after the
        # call so this capability does not change another task's RNG stream.
        random_state = np.random.get_state()
        np.random.seed(parameters.random_seed)
        try:
            scv.pp.normalize_per_cell(adata, counts_per_cell_after=parameters.target_sum,
                                      layers=["spliced", "unspliced"], enforce=True, use_initial_size=False)
            scv.pp.moments(adata, n_neighbors=None, n_pcs=None, mode="connectivities", use_highly_variable=False)
            # In scVelo 0.3.x even False loads an existing HVG mask because
            # the constructor tests "is not None". None avoids that mask.
            scv.tl.velocity(adata, mode="stochastic", filter_genes=False, use_highly_variable=None)
            if "velocity_genes" not in adata.var or not adata.var["velocity_genes"].any():
                raise ValueError("The stochastic kinetic model found no usable velocity genes")
            scv.tl.velocity_graph(adata, n_jobs=1, show_progress_bar=False)
            scv.tl.terminal_states(adata, groupby=None, random_state=parameters.random_seed)
            scv.tl.velocity_pseudotime(adata, groupby=None,
                                       root_key="root_cells", end_key="end_points",
                                       n_dcs=min(parameters.n_dcs, data.n_obs - 2))
            scv.tl.velocity_confidence(adata)
        finally:
            np.random.set_state(random_state)
        adata.X = original_x
        for name in ("n_counts", "n_genes", "n_counts_spliced", "n_counts_unspliced"):
            if name in data.obs:
                adata.obs[name] = data.obs[name].to_numpy()
        if adata.shape != data.shape:
            raise ValueError("RNA velocity must preserve all input cells and genes")
        names = ["velocity_pseudotime", "root_cells", "end_points", "velocity_confidence"]
        for name in names:
            if name not in adata.obs or not np.isfinite(adata.obs[name].to_numpy(dtype=float)).all():
                raise ValueError(f"Stochastic RNA velocity returned missing or nonfinite {name}")
        ids = data.obs["cell_id"].astype(str).to_numpy() if "cell_id" in data.obs else data.obs.index.astype(str).to_numpy()
        if len(set(ids)) != data.n_obs:
            raise ValueError("RNA velocity requires unique cell identifiers")
        table = pd.DataFrame({
            "cell_id": ids, "cluster": data.obs["cluster"].astype(str).to_numpy(),
            "velocity_pseudotime": adata.obs["velocity_pseudotime"].to_numpy(dtype=float),
            "root_score": adata.obs["root_cells"].to_numpy(dtype=float),
            "end_score": adata.obs["end_points"].to_numpy(dtype=float),
            "velocity_confidence": adata.obs["velocity_confidence"].to_numpy(dtype=float),
            "analysis_scope": "within_dataset_model_inferred", "direction_source": "stochastic_kinetic_model",
            "time_units": "unitless_model_pseudotime",
        })
        metrics = {
            "n_cells": data.n_obs, "n_genes": data.n_vars, "n_velocity_genes": int(adata.var["velocity_genes"].sum()),
            "mode": "stochastic", "direction_source": "stochastic_kinetic_model", "root_source": "kinetic_terminal_states",
            "analysis_scope": "within_dataset_model_inferred", "counts_preserved": True,
            "feature_filtering": False, "n_jobs": 1, "summary_groupby": "cluster",
            "kinetic_groupby": "none_global_graph", "confirmatory_allowed": False,
            "mean_velocity_confidence": float(table["velocity_confidence"].mean()),
            "pseudotime_min": float(table["velocity_pseudotime"].min()), "pseudotime_max": float(table["velocity_pseudotime"].max()),
        }
        adata.uns["rna_velocity"] = {
            "method": self.implementation_id, **metrics, "model_assumptions": list(MODEL_ASSUMPTIONS),
            "counts_layers": {"spliced": "spliced_counts", "unspliced": "unspliced_counts", "expression": "counts"},
            "neighbors_source": "native_project_clustering", "caller_root_used": False,
            "n_dcs": min(parameters.n_dcs, data.n_obs - 2), "target_sum": parameters.target_sum,
        }
        sid = ArtifactURI.parse(in_uri).study_id
        outputs = [f"adata://{sid}/rna_velocity/v1", f"table://{sid}/rna_velocity_summary/v1"]
        if contract.expected_outputs and contract.expected_outputs != outputs:
            raise ValueError("RNA velocity requires its canonical paired AnnData/summary outputs")
        versions = {"scvelo": getattr(scv, "__version__", "unknown")}
        for package in ("anndata", "scanpy", "numpy", "scipy"):
            versions[package] = metadata.version(package)
        registry.register(outputs[0], SCData.from_anndata(adata), ArtifactType.ANNDATA, sid, contract.task_id,
                          self.implementation_id, parent_uris=[in_uri], parameters=parameters.model_dump(),
                          software_versions=versions, random_seed=parameters.random_seed, summary_metrics=metrics)
        registry.register(outputs[1], table, ArtifactType.TABLE, sid, contract.task_id,
                          "summarize_stochastic_rna_velocity", parent_uris=[outputs[0]], parameters=parameters.model_dump(),
                          software_versions=versions, random_seed=parameters.random_seed, summary_metrics=metrics)
        return TaskResult(
            task_id=contract.task_id, capability=self.capability_name, method_used=self.implementation_id,
            status=TaskStatus.SUCCESS, input_artifacts=[in_uri], output_artifacts=outputs,
            executed_operations=list(self.contract_operations), metrics=metrics,
        )
