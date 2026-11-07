"""Exercise the real velocity extra on count-valued splicing kinetics.

The model never receives simulated times, kinetic parameters, cell types, or
chosen roots. This integration check verifies the actual optional API and
audited matrix/summary provenance; it does not certify general fate accuracy.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

pytest.importorskip("scvelo")
pytest.importorskip("scanpy")
pytest.importorskip("leidenalg")

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor import ScientificAuditor
from eacbp.capabilities import create_default_capability_registry
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskStatus


def _kinetic_counts():
    rng = np.random.default_rng(31)
    time = np.linspace(0.1, 12, 180)
    spliced, unspliced = [], []
    for gene in range(32):
        alpha, beta, gamma = 25 + gene * 1.5, 1.0, 0.25 + gene % 6 * 0.05
        switch = 5 + gene % 8 * 0.55
        tau = np.minimum(time, switch)
        u = alpha / beta * (1 - np.exp(-beta * tau))
        s = alpha / gamma * (1 - np.exp(-gamma * tau)) + alpha / (gamma - beta) * (np.exp(-gamma * tau) - np.exp(-beta * tau))
        after = time > switch
        delta = time[after] - switch
        u0, s0 = u[after], s[after]
        u[after] = u0 * np.exp(-beta * delta)
        s[after] = s0 * np.exp(-gamma * delta) - u0 * beta / (gamma - beta) * (np.exp(-gamma * delta) - np.exp(-beta * delta))
        spliced.append(rng.poisson(np.maximum(s, 0) + 12))
        unspliced.append(rng.poisson(np.maximum(u, 0) + 4))
    # Stable background genes keep library-size changes from becoming the
    # only signal. These are genuine random count observations in the fixture.
    for gene in range(16):
        spliced.append(rng.poisson(50 + gene, len(time)))
        unspliced.append(rng.poisson(12 + gene % 4, len(time)))
    s = sparse.csr_matrix(np.column_stack(spliced), dtype=np.float32)
    u = sparse.csr_matrix(np.column_stack(unspliced), dtype=np.float32)
    ids = [f"cell{i:03d}" for i in range(len(time))]
    genes = [f"gene{i:02d}" for i in range(s.shape[1])]
    return SCData(
        s, pd.DataFrame({"cell_id": ids}, index=ids), pd.DataFrame({"gene_name": genes}, index=genes),
        layers={"counts": s, "spliced": s, "unspliced": u},
        uns={"is_simulated": True, "data_origin": "synthetic_splicing_kinetic_counts"},
    )


def test_real_stochastic_velocity_from_native_neighbors_preserves_counts_and_infers_orientation(tmp_path):
    data = _kinetic_counts()
    registry = ArtifactRegistry(tmp_path / "artifacts")
    raw_uri = "adata://real_velocity/raw/v1"
    registry.register(raw_uri, data, ArtifactType.ANNDATA, "real_velocity", "input", "input")
    capabilities = create_default_capability_registry()
    normalized = capabilities.execute_contract(TaskContract(
        task_id="normalization", capability="normalization", input_artifacts=[raw_uri],
        parameters={"n_top_genes": data.n_vars},
    ), registry)
    assert normalized.status == TaskStatus.SUCCESS, normalized
    clustered = capabilities.execute_contract(TaskContract(
        task_id="clustering", capability="clustering", method="scanpy_leiden_umap_v1",
        input_artifacts=normalized.output_artifacts, parameters={"random_seed": 13, "n_neighbors": 18, "n_components": 8},
    ), registry)
    assert clustered.status == TaskStatus.SUCCESS, clustered
    task = TaskContract(task_id="velocity", capability="rna_velocity", input_artifacts=clustered.output_artifacts,
                        parameters={"random_seed": 13, "n_dcs": 8})
    capabilities.prepare_contract(task, registry)
    # Keep actual library tracebacks visible in this integration gate.
    result = capabilities.get("rna_velocity").execute(task, registry)
    assert result.status == TaskStatus.SUCCESS, result
    output = registry.get_scdata(result.output_artifacts[0])[1]
    assert output.shape == data.shape
    assert output.uns["velocity_params"]["mode"] == "stochastic"
    assert result.metrics["n_velocity_genes"] >= 2
    assert np.any(np.abs(output.layers["velocity"]) > 0)
    assert output.uns["velocity_graph"].nnz > 0
    assert output.obs["velocity_pseudotime"].max() > output.obs["velocity_pseudotime"].min()
    assert "true_t" not in output.obs and "latent_time" not in output.obs
    report = capabilities.run_validator(task, result, registry)
    assert report.overall_passed, report
    assert ScientificAuditor().audit_task(task, result, registry).overall_passed
    assert "scvelo" in registry.get_metadata(result.output_artifacts[0]).software_versions
    evidence = capabilities.extract_evidence(task, result, report, registry)[0]
    assert evidence.is_simulated is True
    assert evidence.validation_scope == "within_dataset_model_inferred"
    assert evidence.metrics["confirmatory_allowed"] is False
