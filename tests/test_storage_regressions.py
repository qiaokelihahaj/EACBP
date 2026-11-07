"""Regression tests for durable artifact storage and SCData slot fidelity."""

from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.artifact.storage import (
    ArtifactAlreadyExistsError,
    ArtifactDependencyError,
    ArtifactIntegrityError,
    ArtifactSerializationError,
    PayloadSerializer,
)
from eacbp.artifact.uri import ArtifactURI
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType


def _scdata_with_all_slots() -> SCData:
    return SCData(
        X=sparse.csr_matrix(np.arange(12, dtype=np.float32).reshape(3, 4)),
        obs=pd.DataFrame(
            {"cell_id": ["c0", "c1", "c2"], "kind": pd.Categorical(["a", "b", "a"])},
            index=pd.Index(["c0", "c1", "c2"], name="cell"),
        ),
        var=pd.DataFrame({"gene": ["g0", "g1", "g2", "g3"]},
                         index=pd.Index(["g0", "g1", "g2", "g3"], name="feature")),
        layers={"counts": sparse.csc_matrix(np.arange(12, dtype=np.int32).reshape(3, 4))},
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


def _assert_slots_in_fresh_process(
    storage_dir: Path, uri: str, *, block_anndata: bool, uns_values
) -> None:
    block_import = """
import builtins
_original_import = builtins.__import__
def _without_anndata(name, *args, **kwargs):
    if name == 'anndata' or name.startswith('anndata.'):
        raise ImportError('anndata deliberately unavailable for NPZ reload')
    return _original_import(name, *args, **kwargs)
builtins.__import__ = _without_anndata
""" if block_anndata else ""
    script = block_import + f"""
import numpy as np
import pandas as pd
from scipy import sparse
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
data = ArtifactRegistry({json.dumps(str(storage_dir))}).load_scdata({json.dumps(uri)})
assert data.obs.columns.tolist() == ['cell_id', 'kind']
assert data.var.columns.tolist() == ['gene']
assert data.obs.index.tolist() == ['c0', 'c1', 'c2']
assert data.var.index.tolist() == ['g0', 'g1', 'g2', 'g3']
assert data.obs.index.name == 'cell'
assert data.var.index.name == 'feature'
assert data.obs['kind'].astype(str).tolist() == ['a', 'b', 'a']
np.testing.assert_array_equal(data.X.toarray(), np.arange(12, dtype=np.float32).reshape(3, 4))
np.testing.assert_array_equal(data.layers['counts'].toarray(), np.arange(12, dtype=np.int32).reshape(3, 4))
np.testing.assert_array_equal(data.raw['X'], np.arange(15, dtype=np.int32).reshape(3, 5))
assert data.raw['var']['gene'].tolist() == ['r0', 'r1', 'r2', 'r3', 'r4']
np.testing.assert_array_equal(data.obsm['X_pca'], np.arange(6, dtype=np.float32).reshape(3, 2))
np.testing.assert_array_equal(data.obsp['connectivities'].toarray(), np.eye(3))
np.testing.assert_array_equal(data.varm['loadings'], np.arange(8, dtype=np.float32).reshape(4, 2))
np.testing.assert_array_equal(data.varp['gene_graph'], np.eye(4, dtype=np.float32))
np.testing.assert_array_equal(data.uns['nested']['arr'], np.array([1, 2]))
assert list(data.uns['nested']['values']) == {json.dumps(uns_values)}
assert bool(data.uns['flag']) is True
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=Path(__file__).resolve().parents[1],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr or result.stdout


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
    copied.layers["counts"][0, 1] = 99
    assert data.uns["nested"]["arr"][0] == 1
    assert data.layers["counts"][0, 1] == 1

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
    restored = registry.load_scdata(metadata.uri)
    assert restored.uns["nested"]["arr"].tolist() == [1, 2]
    assert np.array_equal(restored.layers["counts"].toarray(), data.layers["counts"].toarray())
    assert restored.raw["X"].shape == (3, 5)
    assert restored.obsp["connectivities"].shape == (3, 3)


@pytest.mark.parametrize("as_sparse", [False, True])
@pytest.mark.parametrize("values", [
    np.array([[2 ** 24 + 1, 3], [5, 2 ** 53 + 1]], dtype=np.int64),
    np.array([[np.nextafter(1.0, 2.0), 0.1], [0.0, np.nextafter(2.0, 3.0)]], dtype=np.float64),
], ids=["integer_counts", "fractional_expression"])
def test_scdata_preserves_x_precision_through_copies_subsets_and_npz(tmp_path, monkeypatch, as_sparse, values):
    def no_anndata():
        raise ArtifactDependencyError("anndata disabled for precision regression")

    monkeypatch.setattr(PayloadSerializer, "_get_anndata_module", no_anndata)
    matrix = sparse.csr_matrix(values) if as_sparse else values.copy()
    data = SCData(matrix, pd.DataFrame(index=["c1", "c2"]), pd.DataFrame(index=["g1", "g2"]))
    # The stored value must not alias the caller's mutable input.
    matrix[0, 0] = 0
    registry = ArtifactRegistry(tmp_path)
    metadata = registry.register("adata://precision/raw/v1", data, ArtifactType.ANNDATA,
                                 "precision", "ingest", "import")
    assert metadata.storage_path.endswith(".npz")
    for output, expected in (
        (data, values),
        (data.copy(), values),
        (SCData.from_dict(data.to_dict()), values),
        (data.subset_obs(np.array([True, False])), values[:1]),
        (data.subset_var(np.array([True, False])), values[:, :1]),
        (ArtifactRegistry(tmp_path).load_scdata(metadata.uri), values),
    ):
        actual = output.X.toarray() if sparse.issparse(output.X) else output.X
        assert actual.dtype == values.dtype
        np.testing.assert_array_equal(actual, expected)


def test_npz_scdata_reloads_in_another_process_without_anndata(tmp_path, monkeypatch):
    def no_anndata():
        raise ArtifactDependencyError("anndata disabled for NPZ storage regression")

    monkeypatch.setattr(PayloadSerializer, "_get_anndata_module", no_anndata)
    registry = ArtifactRegistry(str(tmp_path))
    metadata = registry.register(
        "adata://study/npz_slots/v1", _scdata_with_all_slots(), ArtifactType.ANNDATA,
        "study", "task", "write",
    )
    assert metadata.storage_path.endswith(".npz")

    # Generic reads keep their historical dictionary representation, while
    # the typed boundary normalizes the same verified bytes to SCData.
    generic = registry.load_payload(metadata.uri)
    assert isinstance(generic, dict)
    assert isinstance(registry.get(metadata.uri)[1], dict)
    assert isinstance(registry.load_scdata(metadata.uri), SCData)
    _assert_slots_in_fresh_process(
        tmp_path, metadata.uri, block_anndata=True, uns_values=["a", 2]
    )


def test_h5ad_scdata_reloads_all_slots_in_another_process(tmp_path):
    pytest.importorskip("anndata")
    registry = ArtifactRegistry(str(tmp_path))
    source = _scdata_with_all_slots()
    source.uns["nested"]["values"] = ["a", "b"]
    metadata = registry.register(
        "adata://study/h5ad_slots/v1", source, ArtifactType.ANNDATA,
        "study", "task", "write",
    )
    assert metadata.storage_path.endswith(".h5ad")
    _assert_slots_in_fresh_process(
        tmp_path, metadata.uri, block_anndata=False, uns_values=["a", "b"]
    )


def test_scdata_read_rejects_wrong_type_malformed_payload_and_tampering(tmp_path, monkeypatch):
    registry = ArtifactRegistry(str(tmp_path))
    table = registry.register(
        "table://study/not_scdata/v1", pd.DataFrame({"value": [1]}), ArtifactType.TABLE,
        "study", "task", "write",
    )
    with pytest.raises(TypeError, match="not single-cell data"):
        registry.load_scdata(table.uri)

    malformed = registry.register(
        "adata://study/malformed/v1", {"X": np.ones((2, 2)), "obs": pd.DataFrame(index=["a", "b"])},
        ArtifactType.ANNDATA, "study", "task", "write",
    )
    with pytest.raises(ValueError, match="missing fields.*var"):
        registry.load_scdata(malformed.uri)

    def no_anndata():
        raise ArtifactDependencyError("anndata disabled for NPZ storage regression")

    monkeypatch.setattr(PayloadSerializer, "_get_anndata_module", no_anndata)
    tampered = registry.register(
        "adata://study/tampered/v1", _scdata_with_all_slots(), ArtifactType.ANNDATA,
        "study", "task", "write",
    )
    Path(tampered.storage_path).write_bytes(Path(tampered.storage_path).read_bytes() + b"tamper")
    with pytest.raises(ArtifactIntegrityError):
        registry.load_scdata(tampered.uri)


def test_transaction_typed_reads_cover_parent_and_local_artifacts(tmp_path, monkeypatch):
    def no_anndata():
        raise ArtifactDependencyError("anndata disabled for NPZ transaction regression")

    monkeypatch.setattr(PayloadSerializer, "_get_anndata_module", no_anndata)
    parent = ArtifactRegistry(str(tmp_path))
    source = parent.register(
        "adata://study/source/v1", _scdata_with_all_slots(), ArtifactType.ANNDATA,
        "study", "input", "write",
    )
    transaction = TaskArtifactTransaction(parent)
    try:
        parent_meta, parent_data = transaction.get_scdata(source.uri)
        assert parent_meta.uri == source.uri
        assert parent_data.obs.index.name == "cell"
        assert isinstance(transaction.load_scdata(source.uri), SCData)

        local = transaction.register(
            "adata://study/staged/v1", _scdata_with_all_slots(), ArtifactType.ANNDATA,
            "study", "task", "write",
        )
        local_meta, local_data = transaction.get_scdata(local.uri)
        assert local_meta.uri == local.uri
        assert local_data.shape == (3, 4)
        assert isinstance(transaction.load_scdata(local.uri), SCData)
    finally:
        transaction.close()
