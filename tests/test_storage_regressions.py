"""Regression tests for durable artifact storage and SCData slot fidelity."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.storage import (
    ArtifactAlreadyExistsError,
    ArtifactIntegrityError,
    ArtifactSerializationError,
)
from eacbp.artifact.uri import ArtifactURI
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType


def _scdata_with_all_slots() -> SCData:
    return SCData(
        X=np.arange(12, dtype=np.float32).reshape(3, 4),
        obs=pd.DataFrame(
            {"cell_id": ["c0", "c1", "c2"], "kind": pd.Categorical(["a", "b", "a"])}
        ),
        var=pd.DataFrame({"gene": ["g0", "g1", "g2", "g3"]}),
        layers={"counts": np.arange(12, dtype=np.int32).reshape(3, 4)},
        raw={
            "X": np.arange(15, dtype=np.int32).reshape(3, 5),
            "var": pd.DataFrame({"gene": ["r0", "r1", "r2", "r3", "r4"]}),
        },
        obsm={"X_pca": np.arange(6, dtype=np.float32).reshape(3, 2)},
        obsp={"connectivities": sparse.eye(3, format="csr")},
        varm={"loadings": np.arange(8, dtype=np.float32).reshape(4, 2)},
        varp={"gene_graph": np.eye(4, dtype=np.float32)},
        uns={"nested": {"arr": np.array([1, 2]), "values": ["a", 2]}, "flag": True},
    )


def test_registry_persists_metadata_lineage_and_provenance(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    parent = registry.register(
        "adata://study/raw/v1",
        _scdata_with_all_slots().to_dict(),
        ArtifactType.ANNDATA,
        "study",
        "ingest",
        "ingest",
        summary_metrics={"is_simulated": True, "data_origin": "synthetic"},
    )
    child = registry.register(
        "table://study/summary/v1",
        pd.DataFrame({"value": [1]}),
        ArtifactType.TABLE,
        "study",
        "summary",
        "summarize",
        parent_uris=[parent.uri],
    )

    rebuilt = ArtifactRegistry(str(tmp_path))
    loaded_meta, loaded_table = rebuilt.get(child.uri)
    assert loaded_table.loc[0, "value"] == 1
    assert loaded_meta.summary_metrics["is_simulated"] is True
    assert loaded_meta.summary_metrics["data_origin"] == "synthetic"
    assert parent.uri in rebuilt.lineage.get_ancestors(child.uri)
    assert rebuilt.exists(parent.uri)


def test_read_rejects_tampered_payload(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    metadata = registry.register(
        "json://study/report/v1",
        {"answer": 42},
        ArtifactType.JSON,
        "study",
        "task",
        "write",
    )
    Path(metadata.storage_path).write_text('{"answer": 43}', encoding="utf-8")
    with pytest.raises(ArtifactIntegrityError):
        registry.load_payload(metadata.uri)


def test_failed_structured_write_leaves_no_target_and_no_pickle(tmp_path):
    class Unsupported:
        pass

    registry = ArtifactRegistry(str(tmp_path))
    with pytest.raises(ArtifactSerializationError):
        registry.register(
            "adata://study/bad/v1",
            {"X": np.ones((1, 1)), "unsupported": Unsupported()},
            ArtifactType.ANNDATA,
            "study",
            "task",
            "write",
        )
    assert not list((tmp_path / "study" / "bad").glob("v1.*"))


def test_cross_instance_registration_has_one_winner(tmp_path):
    def register_one(value: int):
        registry = ArtifactRegistry(str(tmp_path))
        try:
            return registry.register(
                "json://study/concurrent/v1",
                {"value": value},
                ArtifactType.JSON,
                "study",
                f"task-{value}",
                "write",
            ).summary_metrics
        except ArtifactAlreadyExistsError:
            return "exists"

    with ThreadPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(register_one, range(6)))
    assert results.count("exists") == 5
    assert ArtifactRegistry(str(tmp_path)).get("json://study/concurrent/v1")[1]["value"] in range(6)


def test_uri_and_lineage_boundaries(tmp_path):
    with pytest.raises(ValueError):
        ArtifactURI.parse("json://study/../escape/v1")
    registry = ArtifactRegistry(str(tmp_path))
    with pytest.raises(ValueError, match="itself"):
        registry.register(
            "json://study/a/v1", {}, ArtifactType.JSON, "study", "task", "write", ["json://study/a/v1"]
        )
    registry.register("json://study/a/v1", {}, ArtifactType.JSON, "study", "task", "write", ["json://study/b/v1"])
    with pytest.raises(ValueError, match="cycle"):
        registry.register("json://study/b/v1", {}, ArtifactType.JSON, "study", "task", "write", ["json://study/a/v1"])


def test_scdata_copy_subset_and_safe_fallback_preserve_slots(tmp_path):
    data = _scdata_with_all_slots()
    copied = data.copy()
    copied.uns["nested"]["arr"][0] = 99
    copied.layers["counts"][0, 0] = 99
    assert data.uns["nested"]["arr"][0] == 1
    assert data.layers["counts"][0, 0] == 0

    subset = data.subset_obs(np.array([True, False, True]))
    assert subset.shape == (2, 4)
    assert subset.layers["counts"].shape == (2, 4)
    assert subset.obsp["connectivities"].shape == (2, 2)
    assert subset.raw["X"].shape == (2, 5)

    variable_subset = data.subset_var(np.array([True, False, True, False]))
    assert variable_subset.shape == (3, 2)
    assert variable_subset.varm["loadings"].shape == (2, 2)
    assert variable_subset.varp["gene_graph"].shape == (2, 2)
    assert variable_subset.raw["X"].shape == (3, 5)

    registry = ArtifactRegistry(str(tmp_path))
    metadata = registry.register(
        "adata://study/slots/v1",
        data.to_dict(),
        ArtifactType.ANNDATA,
        "study",
        "task",
        "write",
    )
    assert metadata.storage_path.endswith(".npz") or metadata.storage_path.endswith(".h5ad")
    restored = registry.load_payload(metadata.uri)
    restored = restored if isinstance(restored, SCData) else SCData.from_dict(restored)
    assert restored.uns["nested"]["arr"].tolist() == [1, 2]
    assert np.array_equal(restored.layers["counts"], data.layers["counts"])
    assert restored.raw["X"].shape == (3, 5)
    assert restored.obsp["connectivities"].shape == (3, 3)
