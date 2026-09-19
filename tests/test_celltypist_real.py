"""A synthetic model verifies library integration, not biological annotation accuracy."""
import numpy as np
import pandas as pd
import pytest
pytest.importorskip("celltypist")


def test_actual_celltypist_local_model_roundtrip(tmp_path, monkeypatch):
    import celltypist
    from eacbp.capabilities.advanced_qc import CellTypistAnnotationCapability
    from eacbp.capabilities.sc_data import SCData
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.schemas.artifact import ArtifactType
    from eacbp.schemas.task import TaskContract
    rng = np.random.default_rng(6)
    counts = rng.poisson(3, (60, 30))
    counts[:30, :10] += 10
    counts[30:, 10:20] += 10
    x = np.log1p(counts / counts.sum(axis=1)[:, None] * 10000)
    genes = [f"g{i}" for i in range(30)]
    labels = ["SyntheticA"] * 30 + ["SyntheticB"] * 30
    model = celltypist.train(X=x, labels=labels, genes=genes, n_jobs=1, max_iter=100, check_expression=True, use_SGD=True)
    path = tmp_path / "synthetic_model.pkl"
    model.write(str(path))
    def reject_remote_discovery(*args, **kwargs):
        raise AssertionError("Local annotation must not discover/download remote models")
    monkeypatch.setattr(celltypist.models, "get_all_models", reject_remote_discovery)
    data = SCData(x, pd.DataFrame({"cell_type": ["ExistingLabel"] * 60}), pd.DataFrame(index=genes), layers={"counts": counts})
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "adata://toy/raw/v1"
    registry.register(uri, data, ArtifactType.ANNDATA, "toy", "input", "input", summary_metrics={"is_simulated": True})
    contract = TaskContract(task_id="annotate", capability="cell_annotation", input_artifacts=[uri],
                            expected_outputs=["adata://toy/annotated/v1"], parameters={"model_path": str(path), "model_species": "synthetic", "model_tissue": "synthetic"})
    result = CellTypistAnnotationCapability().execute(contract, registry)
    _, output = registry.get(result.output_artifacts[0])
    if isinstance(output, dict):
        output = SCData.from_dict(output)
    assert output.obs.cell_type.eq("ExistingLabel").all()
    assert set(output.obs.cell_type_celltypist) <= {"SyntheticA", "SyntheticB"}
    assert output.obs.cell_type_celltypist_conflict.all()
    np.testing.assert_array_equal(output.layers["counts"], counts)
    from eacbp.auditor import ScientificAuditor
    contract.validation_requirements = ["advanced_qc_integrity"]
    audited = ScientificAuditor().audit_task(contract, result, registry)
    assert audited.overall_passed, [(c.check_name, c.message) for c in audited.checks if not c.passed]
