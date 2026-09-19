import numpy as np
import pandas as pd
import pytest
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.perturbation.compound import CompoundPerturbationCapability
from eacbp.capabilities.perturbation.genetic import GeneticPerturbationCapability
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract


def test_deg_table_does_not_manufacture_counterfactual_cells(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    uri = "table://p/deg/v1"
    reg.register(uri, pd.DataFrame({"gene": ["G"]}), ArtifactType.TABLE, "p", "input", "input")
    task = TaskContract(task_id="c", capability="compound_perturbation_simulation", input_artifacts=[uri])
    with pytest.raises(ValueError):
        CompoundPerturbationCapability().execute(task, reg)
    assert len(reg.list_artifacts()) == 1


def test_single_cohort_cannot_define_disease_signature(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    data = SCData(np.ones((8, 4)), pd.DataFrame({"condition": ["one"] * 8}), pd.DataFrame({"gene_name": ["Apoe", "Trem2", "P2ry12", "G"]}))
    uri = "adata://p/raw/v1"
    reg.register(uri, data, ArtifactType.ANNDATA, "p", "input", "input")
    task = TaskContract(task_id="c", capability="compound_perturbation_simulation", input_artifacts=[uri], parameters={"compound_name": "Bexarotene"})
    with pytest.raises(ValueError, match="Disease signature requires"):
        CompoundPerturbationCapability().execute(task, reg)


def test_genetic_projection_does_not_mix_unrelated_embedding_basis(tmp_path):
    x = np.random.default_rng(3).poisson(3, (20, 4)).astype(float)
    outputs = []
    for i, embedding in enumerate([np.ones((20, 2)), np.ones((20, 2)) * 999]):
        reg = ArtifactRegistry(str(tmp_path / str(i)))
        data = SCData(x, pd.DataFrame({"condition": ["control"] * 10 + ["AD"] * 10}),
                      pd.DataFrame({"gene_name": ["Apoe", "Trem2", "G1", "G2"]}), obsm={"X_pca": embedding})
        uri = "adata://p/raw/v1"
        reg.register(uri, data, ArtifactType.ANNDATA, "p", "input", "input")
        task = TaskContract(task_id="g", capability="genetic_perturbation_simulation", input_artifacts=[uri], parameters={"target_gene": "Trem2"})
        result = GeneticPerturbationCapability().execute(task, reg)
        outputs.append(reg.get(result.output_artifacts[0])[1].obsm["X_pca"])
    assert np.allclose(outputs[0], outputs[1])
