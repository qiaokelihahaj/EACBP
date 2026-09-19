"""Contract-level tests for real optional QC/annotation adapters.

The tests that exercise the adapter wiring use tiny in-process doubles for
optional third-party libraries.  The real Scrublet test is marked optional and
is run by environments that install Scanpy's scrublet extra.
"""

from __future__ import annotations

import hashlib
import sys
import types

import numpy as np
import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.advanced_qc import (
    CellBenderBackgroundRemovalCapability,
    CellTypistAnnotationCapability,
    ScrubletDoubletCapability,
    _map_external_matrix_to_target,
)
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskStatus


def _register_input(tmp_path, *, n_cells=8, n_genes=5, existing_types=True):
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    obs = pd.DataFrame(
        {
            "cell_id": [f"cell_{i}" for i in range(n_cells)],
            "batch": ["library_a" if i < n_cells // 2 else "library_b" for i in range(n_cells)],
        }
    )
    if existing_types:
        obs["cell_type"] = ["Existing"] * n_cells
    var = pd.DataFrame({"gene_name": [f"gene_{i}" for i in range(n_genes)]})
    counts = np.arange(n_cells * n_genes, dtype=np.int64).reshape(n_cells, n_genes)
    counts[0, 0] = 1
    data = SCData(counts, obs, var, layers={"counts": counts.copy()})
    uri = "adata://advanced/raw/v1"
    registry.register(uri, data, ArtifactType.ANNDATA, "advanced", "input", "input")
    return registry, uri, data


def test_capability_interfaces_are_explicit():
    scrublet = ScrubletDoubletCapability()
    celltypist = CellTypistAnnotationCapability()
    cellbender = CellBenderBackgroundRemovalCapability()
    assert (scrublet.capability_name, scrublet.implementation_id) == (
        "doublet_detection",
        "scanpy_scrublet_v1",
    )
    assert (celltypist.capability_name, celltypist.implementation_id) == (
        "cell_annotation",
        "celltypist_local_v1",
    )
    assert (cellbender.capability_name, cellbender.implementation_id) == (
        "background_removal",
        "cellbender_cli_v1",
    )


def test_scrublet_marks_by_batch_preserves_counts_and_honors_expected_output(tmp_path, monkeypatch):
    registry, uri, source = _register_input(tmp_path)

    def fake_scrublet(adata, **kwargs):
        adata.obs["doublet_score"] = np.linspace(0.1, 0.8, adata.n_obs)
        adata.obs["predicted_doublet"] = np.array([False, True, False, False, True, False, False, False])
        adata.uns["scrublet"] = {
            "batches": {
                "library_a": {"threshold": 0.3},
                "library_b": {"threshold": 0.4},
            },
            "parameters": {"expected_doublet_rate": kwargs["expected_doublet_rate"]},
        }

    fake_scanpy = types.ModuleType("scanpy")
    fake_scanpy.pp = types.SimpleNamespace(scrublet=fake_scrublet)
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setitem(sys.modules, "scrublet", types.ModuleType("scrublet"))

    task = TaskContract(
        task_id="doublet",
        capability="doublet_detection",
        input_artifacts=[uri],
        expected_outputs=["adata://advanced/doublets/v7", "table://advanced/doublet_scores/v1"],
        parameters={"n_prin_comps": 4, "filter_doublets": False},
    )
    result = ScrubletDoubletCapability().execute(task, registry)
    assert result.status == TaskStatus.SUCCESS
    assert result.output_artifacts == [
        "adata://advanced/doublets/v7",
        "table://advanced/doublet_scores/v1",
    ]
    _, output = registry.get(result.output_artifacts[0])
    assert output.n_obs == source.n_obs
    assert np.array_equal(output.layers["counts"], source.layers["counts"])
    assert output.obs["predicted_doublet"].sum() == 2
    assert result.metrics["n_cells_before"] == 8
    assert result.metrics["n_cells_after"] == 8
    assert result.metrics["threshold_by_batch"] == {"library_a": 0.3, "library_b": 0.4}
    assert "filter_cells" not in result.executed_operations


def test_scrublet_filter_requires_explicit_parameter_and_records_before_after(tmp_path, monkeypatch):
    registry, uri, source = _register_input(tmp_path)

    def fake_scrublet(adata, **kwargs):
        adata.obs["doublet_score"] = np.array([0.3, 0.1, 0.1, 0.1, 0.4, 0.1, 0.1, 0.1])
        adata.obs["predicted_doublet"] = np.array([True, False, False, False, True, False, False, False])
        adata.uns["scrublet"] = {"batches": {"library_a": {"threshold": 0.2}, "library_b": {"threshold": 0.2}}}

    fake_scanpy = types.ModuleType("scanpy")
    fake_scanpy.pp = types.SimpleNamespace(scrublet=fake_scrublet)
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setitem(sys.modules, "scrublet", types.ModuleType("scrublet"))
    task = TaskContract(
        task_id="doublet_filter",
        capability="doublet_detection",
        input_artifacts=[uri],
        expected_outputs=["adata://advanced/doublets_filtered/v1", "table://advanced/filter_scores/v1"],
        parameters={"filter_doublets": True, "n_prin_comps": 3},
    )
    result = ScrubletDoubletCapability().execute(task, registry)
    _, output = registry.get(result.output_artifacts[0])
    assert output.n_obs == 6
    assert result.metrics["n_cells_before"] == source.n_obs
    assert result.metrics["n_cells_after"] == 6
    assert result.metrics["n_doublets_filtered"] == 2
    assert "filter_cells" in result.executed_operations
    assert output.layers["counts"].shape == (6, source.n_vars)
    from eacbp.auditor.advanced_qc import AdvancedQCValidator
    validator = AdvancedQCValidator()
    assert validator.audit(task, result, registry).overall_passed
    original_get = registry.get
    def tampered_get(artifact_uri):
        metadata, payload = original_get(artifact_uri)
        if artifact_uri == result.output_artifacts[1]:
            payload = payload.copy()
            # A removed cell below the threshold must invalidate the audit.
            payload.loc[0, "doublet_score"] = 0.01
        return metadata, payload
    monkeypatch.setattr(registry, "get", tampered_get)
    assert not validator.audit(task, result, registry).overall_passed


def test_scrublet_real_scanpy_small_counts_by_library_batch(tmp_path):
    """Run the actual Scanpy/Scrublet stack when its optional extra is present."""

    pytest.importorskip("scanpy")
    pytest.importorskip("skimage")
    registry = ArtifactRegistry(str(tmp_path / "real_artifacts"))
    rng = np.random.default_rng(31)
    counts = rng.poisson(2.0, size=(64, 20)).astype(np.int32)
    obs = pd.DataFrame(
        {
            "cell_id": [f"real_cell_{i}" for i in range(64)],
            "batch": ["library_a"] * 32 + ["library_b"] * 32,
        }
    )
    var = pd.DataFrame({"gene_name": [f"real_gene_{i}" for i in range(20)]})
    data = SCData(counts, obs, var, layers={"counts": counts.copy()})
    input_uri = "adata://advanced_real/raw/v1"
    registry.register(input_uri, data, ArtifactType.ANNDATA, "advanced_real", "input", "input")
    task = TaskContract(
        task_id="real_scrublet",
        capability="doublet_detection",
        input_artifacts=[input_uri],
        expected_outputs=["adata://advanced_real/doublets/v1"],
        parameters={
            "batch_key": "batch",
            "threshold": 0.2,
            "use_approx_neighbors": False,
            "n_prin_comps": 10,
            "random_seed": 13,
        },
    )
    result = ScrubletDoubletCapability().execute(task, registry)
    _, output = registry.get(result.output_artifacts[0])
    assert result.method_used == "scanpy_scrublet_v1"
    assert output.n_obs == data.n_obs
    assert {"doublet_score", "predicted_doublet"} <= set(output.obs.columns)
    assert np.array_equal(output.layers["counts"], counts)
    assert result.metrics["n_cells_before"] == 64
    assert result.metrics["n_cells_after"] == 64
    assert set(result.metrics["threshold_by_batch"]) == {"library_a", "library_b"}
    from eacbp.auditor import ScientificAuditor
    task.validation_requirements = ["advanced_qc_integrity"]
    audited = ScientificAuditor().audit_task(task, result, registry)
    assert audited.overall_passed, [(c.check_name, c.message) for c in audited.checks if not c.passed]


def test_scrublet_rejects_fractional_counts_before_optional_runtime(tmp_path):
    registry, uri, _ = _register_input(tmp_path)
    _, payload = registry.get(uri)
    payload.layers["counts"] = payload.layers["counts"].astype(float)
    payload.layers["counts"][0, 0] = 0.5
    # Re-register a separate URI because registry artifacts are immutable.
    registry.register("adata://advanced/fractional/v1", payload, ArtifactType.ANNDATA, "advanced", "input2", "input")
    task = TaskContract(
        task_id="bad_counts",
        capability="doublet_detection",
        input_artifacts=["adata://advanced/fractional/v1"],
        parameters={},
    )
    with pytest.raises(ValueError, match="integer-valued raw counts"):
        ScrubletDoubletCapability().execute(task, registry)


def test_celltypist_local_model_hash_and_existing_labels_are_preserved(tmp_path, monkeypatch):
    registry, uri, source = _register_input(tmp_path, n_cells=3, n_genes=4)
    model_path = tmp_path / "local_model.pkl"
    model_path.write_bytes(b"local-celltypist-model")
    model_hash = hashlib.sha256(model_path.read_bytes()).hexdigest()

    class FakeModel:
        @classmethod
        def load(cls, path):
            assert model_path.__class__(path).resolve() == model_path.resolve()
            return types.SimpleNamespace(metadata={"species": "human", "tissue": "kidney"})

    class FakeAnnotation:
        predicted_labels = pd.DataFrame(
            {"predicted_labels": ["ModelType", "Unknown", "ModelType"], "conf_score": [0.9, 0.2, 0.8]}
        )

    fake_celltypist = types.ModuleType("celltypist")
    fake_celltypist.__path__ = []
    fake_celltypist.annotate = lambda adata, **kwargs: FakeAnnotation()
    fake_models = types.ModuleType("celltypist.models")
    fake_models.Model = FakeModel
    monkeypatch.setitem(sys.modules, "celltypist", fake_celltypist)
    monkeypatch.setitem(sys.modules, "celltypist.models", fake_models)

    task = TaskContract(
        task_id="annotation",
        capability="cell_annotation",
        input_artifacts=[uri],
        expected_outputs=["adata://advanced/annotated/v1", "json://advanced/annotation_report/v1"],
        parameters={
            "model_path": str(model_path),
            "model_metadata": {"species": "human", "tissue": "kidney"},
            "external_resource_sha256": {"model_path": model_hash},
        },
    )
    result = CellTypistAnnotationCapability().execute(task, registry)
    _, output = registry.get(result.output_artifacts[0])
    assert output.obs["cell_type"].tolist() == source.obs["cell_type"].tolist()
    assert output.obs["cell_type_existing"].tolist() == source.obs["cell_type"].tolist()
    assert output.obs["cell_type_celltypist"].tolist() == ["ModelType", "Unknown", "ModelType"]
    assert output.obs["cell_type_celltypist_unknown"].tolist() == [False, True, False]
    assert result.metrics["model_sha256"] == model_hash
    assert result.metrics["model_metadata"]["species"] == "human"


def test_celltypist_rejects_model_hash_mismatch(tmp_path):
    registry, uri, _ = _register_input(tmp_path, n_cells=3, n_genes=4)
    model_path = tmp_path / "local_model.pkl"
    model_path.write_bytes(b"local-celltypist-model")
    task = TaskContract(
        task_id="annotation_bad_hash",
        capability="cell_annotation",
        input_artifacts=[uri],
        parameters={"model_path": str(model_path), "external_resource_sha256": {"model_path": "0" * 64}},
    )
    with pytest.raises(ValueError, match="hash mismatch"):
        CellTypistAnnotationCapability().execute(task, registry)


def test_celltypist_explicit_use_writes_cell_type_without_mutating_cluster(tmp_path, monkeypatch):
    registry, uri, _ = _register_input(tmp_path, n_cells=3, n_genes=4, existing_types=False)
    _, payload = registry.get(uri)
    payload.obs["cluster"] = ["0", "1", "0"]
    registry.register("adata://advanced/clustered/v1", payload, ArtifactType.ANNDATA, "advanced", "cluster", "cluster")
    model_path = tmp_path / "local_model.pkl"
    model_path.write_bytes(b"local-celltypist-model")

    class FakeModel:
        @classmethod
        def load(cls, path):
            return types.SimpleNamespace(metadata={"species": "unknown", "tissue": "unknown"})

    class FakeAnnotation:
        predicted_labels = pd.DataFrame({"predicted_labels": ["A", "B", "A"]})

    fake_celltypist = types.ModuleType("celltypist")
    fake_celltypist.__path__ = []
    fake_celltypist.annotate = lambda adata, **kwargs: FakeAnnotation()
    fake_models = types.ModuleType("celltypist.models")
    fake_models.Model = FakeModel
    monkeypatch.setitem(sys.modules, "celltypist", fake_celltypist)
    monkeypatch.setitem(sys.modules, "celltypist.models", fake_models)
    task = TaskContract(
        task_id="annotation_cluster",
        capability="cell_annotation",
        input_artifacts=["adata://advanced/clustered/v1"],
        expected_outputs=["adata://advanced/cluster_annotation/v1"],
        parameters={
            "model_path": str(model_path),
            "use_as_cell_type": True,
            "model_metadata": {"species": "unknown", "tissue": "unknown"},
        },
    )
    result = CellTypistAnnotationCapability().execute(task, registry)
    _, output = registry.get(result.output_artifacts[0])
    assert output.obs["cluster"].tolist() == ["0", "1", "0"]
    assert output.obs["cell_type"].tolist() == ["A", "B", "A"]


def test_cellbender_requires_explicit_unfiltered_input_executable_and_output(tmp_path):
    registry, uri, _ = _register_input(tmp_path)
    task = TaskContract(
        task_id="cellbender_missing",
        capability="background_removal",
        input_artifacts=[uri],
        parameters={},
    )
    with pytest.raises(ValueError, match="unfiltered_input_path"):
        CellBenderBackgroundRemovalCapability().execute(task, registry)


def test_cellbender_rejects_non_h5_output_path(tmp_path):
    registry, uri, _ = _register_input(tmp_path)
    raw_path = tmp_path / "raw.h5"
    raw_path.write_bytes(b"raw")
    executable = tmp_path / "cellbender.exe"
    executable.write_bytes(b"executable")
    task = TaskContract(
        task_id="cellbender_format",
        capability="background_removal",
        input_artifacts=[uri],
        parameters={
            "unfiltered_input_path": str(raw_path),
            "executable": str(executable),
            "output_path": str(tmp_path / "output.h5ad"),
        },
    )
    with pytest.raises(ValueError, match="must end with .h5"):
        CellBenderBackgroundRemovalCapability().execute(task, registry)


def test_cellbender_mapping_rejects_duplicate_or_missing_target_ids():
    matrix = np.ones((3, 2), dtype=np.float32)
    with pytest.raises(ValueError, match="duplicated"):
        _map_external_matrix_to_target(
            matrix,
            np.array(["cell_a", "cell_a", "cell_b"]),
            np.array(["gene_a", "gene_b"]),
            np.array(["cell_a", "cell_b"]),
            np.array(["gene_a", "gene_b"]),
        )
    with pytest.raises(ValueError, match="missing target cell IDs"):
        _map_external_matrix_to_target(
            matrix,
            np.array(["cell_a", "cell_b", "cell_c"]),
            np.array(["gene_a", "gene_b"]),
            np.array(["cell_a", "cell_missing"]),
            np.array(["gene_a", "gene_b"]),
        )


def test_cellbender_cli_preserves_counts_and_stores_corrected_layer(tmp_path, monkeypatch):
    registry, uri, source = _register_input(tmp_path, n_cells=4, n_genes=3)
    raw_path = tmp_path / "raw_unfiltered.h5"
    raw_path.write_bytes(b"unfiltered droplet input")
    executable = tmp_path / "cellbender.exe"
    executable.write_bytes(b"fake executable")
    output_path = tmp_path / "corrected.h5"

    external_obs = pd.concat(
        [
            pd.DataFrame({"cell_id": ["extra_droplet"], "batch": ["library_extra"]}),
            source.obs.iloc[[3, 1, 0, 2]].copy(),
        ],
        ignore_index=True,
    )
    external_var = pd.concat(
        [
            pd.DataFrame({"gene_name": ["extra_gene"]}),
            source.var.iloc[[2, 0, 1]].copy(),
        ],
        ignore_index=True,
    )
    external_matrix = np.arange(5 * 4, dtype=np.float32).reshape(5, 4)
    external = SCData(external_matrix, external_obs, external_var).to_anndata()
    monkeypatch.setattr(
        CellBenderBackgroundRemovalCapability,
        "_read_output",
        staticmethod(lambda path: external),
    )

    def fake_run(command, **kwargs):
        output_path.write_bytes(b"external corrected output")
        return types.SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr("eacbp.capabilities.advanced_qc.subprocess.run", fake_run)
    task = TaskContract(
        task_id="cellbender",
        capability="background_removal",
        input_artifacts=[uri],
        expected_outputs=["adata://advanced/corrected/v3", "json://advanced/cellbender_report/v1"],
        parameters={
            "unfiltered_input_path": str(raw_path),
            "executable": str(executable),
            "output_path": str(output_path),
            "run_cwd": str(tmp_path),
            "external_resource_sha256": {
                "unfiltered_input_path": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
                "executable": hashlib.sha256(executable.read_bytes()).hexdigest(),
            },
        },
    )
    result = CellBenderBackgroundRemovalCapability().execute(task, registry)
    _, output = registry.get(result.output_artifacts[0])
    assert result.status == TaskStatus.SUCCESS
    assert np.array_equal(output.layers["counts"], source.layers["counts"])
    expected = external_matrix[[3, 2, 4, 1]][:, [2, 3, 1]]
    assert np.array_equal(output.layers["corrected_counts"], expected)
    assert output.uns["cellbender"]["cell_gene_alignment"] == {
        "target_n_cells": 4,
        "target_n_genes": 3,
        "external_n_cells": 5,
        "external_n_genes": 4,
        "mapped_n_cells": 4,
        "mapped_n_genes": 3,
        "extra_external_cells": 1,
        "extra_external_genes": 1,
        "cell_order_exact": False,
        "gene_order_exact": False,
        "cell_ids_unique": True,
        "gene_ids_unique": True,
    }
    assert result.metrics["run_cwd"] == str(tmp_path.resolve())
    # This is an adapter-boundary test: the fake CLI output is intentionally
    # not claimed to be a real CellBender matrix.  The independent auditor must
    # reject it and, after mutation, report the changed external hash.
    from eacbp.auditor.advanced_qc import AdvancedQCValidator

    output_path.write_bytes(b"mutated external output")
    audit = AdvancedQCValidator().audit(task, result, registry)
    assert not audit.overall_passed
    assert any("output_path hash differs" in check.message for check in audit.checks if not check.passed)
