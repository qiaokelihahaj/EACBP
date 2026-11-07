"""New Table v1 binary semantics and unchanged legacy CSV compatibility."""
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from eacbp.artifact.portability import export_artifact_bundle, import_artifact_bundle
from eacbp.artifact.registry import ArtifactRegistry, ArtifactAuditAccessError
from eacbp.artifact.storage import (
    ArtifactAlreadyExistsError, ArtifactSerializationError, ArtifactStorageBackend,
    PayloadSerializer,
)
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.auditor.base import ValidationReport
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def sample_table():
    # P4 stored token 12538.393112853733 was loaded by the default CSV parser
    # as 12538.393112853732. It is 36775107 / 2933, exactly in float64.
    return pd.DataFrame({
        "mean": np.array([36775107 / 2933, np.nan, -0.0, np.inf], dtype=np.float64),
        "float32": np.array([np.float32(1 / 7), np.nan, -0.0, np.inf], dtype=np.float32),
        "large_int": np.array([2**53 + 1, np.iinfo(np.int64).max, -1, 0], dtype=np.int64),
        "large_uint": np.array([np.iinfo(np.uint64).max, 2**53 + 1, 1, 0], dtype=np.uint64),
        "text": pd.Series(["001", "NA", "", "NULL"], dtype=object),
        "object_missing": pd.Series([None, np.nan, pd.NA, "nan"], dtype=object),
        "nullable_integer": pd.Series([2**53 + 1, pd.NA, -1, 0], dtype="Int64"),
        "nullable_float": pd.Series([1.0, pd.NA, -0.0, 1 / 7], dtype="Float32"),
        "nullable_bool": pd.Series([True, pd.NA, False, True], dtype="boolean"),
        "nullable_text": pd.Series(["001", pd.NA, "", "NA"], dtype="string[python]"),
        "category": pd.Categorical(["NA", None, "001", ""], categories=["", "001", "NA"], ordered=True),
    }).set_axis(pd.Index(["001", "NA", "", "NULL"], dtype=object, name="cell_id"))


def assert_exact_table(expected, observed):
    pd.testing.assert_frame_equal(expected, observed, check_exact=True)
    for column in ("mean", "float32"):
        if column in expected:
            assert expected[column].to_numpy().tobytes() == observed[column].to_numpy().tobytes()
    if "object_missing" in expected:
        assert observed.object_missing.iloc[0] is None
        assert type(observed.object_missing.iloc[1]) is float and np.isnan(observed.object_missing.iloc[1])
        assert observed.object_missing.iloc[2] is pd.NA
        assert observed.object_missing.iloc[3] == "nan"


def save_table(root, frame=None):
    registry = ArtifactRegistry(str(root))
    frame = sample_table() if frame is None else frame
    metadata = registry.register("table://s/table/v1", frame, ArtifactType.TABLE, "s", "task", "compute")
    return registry, metadata


def read_contract(path):
    with np.load(path, allow_pickle=False) as stored:
        return json.loads(stored[PayloadSerializer.TABLE_CONTRACT_KEY].tobytes().decode("utf-8"))


def test_binary_table_preserves_values_dtypes_missing_and_axes(tmp_path):
    frame = sample_table()
    registry, metadata = save_table(tmp_path, frame)
    assert metadata.storage_path.endswith(".table.npz")
    _, loaded = ArtifactRegistry(str(tmp_path)).get(metadata.uri)
    assert_exact_table(frame, loaded)
    assert loaded.loc["001", "mean"].hex() == "0x1.87d325185a136p+13"
    contract = read_contract(metadata.storage_path)
    assert (contract["format"], contract["version"], contract["decoder"]) == ("eacbp.table.npz", 1, "eacbp.table.v1")
    assert contract["shape"] == [4, 11]
    assert contract["column_values"][0]["dtype"] == np.dtype("float64").str
    assert contract["column_values"][1]["dtype"] == np.dtype("float32").str
    assert registry.storage.exists(metadata.uri, metadata.type)
    assert_exact_table(frame, registry.storage.load(metadata.uri, metadata.type, metadata.sha256_hash))


@pytest.mark.parametrize("dtype", ["int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64", "float16", "float32", "float64", "bool",
                                    "Int8", "Int16", "Int32", "Int64", "UInt8", "UInt16", "UInt32", "UInt64", "Float32", "Float64", "boolean"])
def test_native_and_nullable_numeric_dtypes_are_preserved(tmp_path, dtype):
    nullable = dtype[:1].isupper() or dtype == "boolean"
    expected = pd.DataFrame({"value": pd.Series([1, pd.NA if nullable else 0, 1], dtype=dtype)})
    registry, metadata = save_table(tmp_path, expected)
    pd.testing.assert_frame_equal(expected, registry.get(metadata.uri)[1], check_exact=True)


def test_nullable_float_distinguishes_unmasked_nan_from_masked_na(tmp_path):
    values = pd.arrays.FloatingArray(np.array([np.nan, -0.0, 1.0]), np.array([False, True, False]))
    expected = pd.DataFrame({"value": values})
    registry, metadata = save_table(tmp_path, expected)
    actual = registry.get(metadata.uri)[1]
    pd.testing.assert_frame_equal(expected, actual, check_exact=True)
    assert actual.value.isna().tolist() == [False, True, False]
    assert np.isnan(actual.value.iloc[0]) and actual.value.iloc[1] is pd.NA


def test_object_numpy_scalar_types_and_text_identity_are_preserved(tmp_path):
    values = [np.str_("001"), np.str_("NA"), np.float32(1 / 7), np.int16(7), np.bool_(True)]
    expected = pd.DataFrame({"value": pd.Series(values, dtype=object)})
    registry, metadata = save_table(tmp_path, expected)
    actual = registry.get(metadata.uri)[1]
    pd.testing.assert_frame_equal(expected, actual, check_exact=True)
    assert [type(value) for value in actual.value] == [type(value) for value in values]


@pytest.mark.parametrize("missing", [pd.NA, np.nan])
def test_string_dtype_missing_policy_is_explicit(tmp_path, missing):
    expected = pd.DataFrame({"id": pd.Series(["001", "NA", "", missing], dtype=pd.StringDtype(storage="python", na_value=missing))})
    registry, metadata = save_table(tmp_path, expected)
    pd.testing.assert_frame_equal(expected, registry.get(metadata.uri)[1], check_exact=True)


@pytest.mark.parametrize("index", [pd.RangeIndex(2, 14, 3, name="rows"), pd.Index([1, 1, 2, 3], name="rows"),
                                  pd.CategoricalIndex(["001", "NA", "", "001"], categories=["NA", "001", ""], ordered=True, name="rows"),
                                  pd.Index(pd.array(["001", pd.NA, "", "NA"], dtype="string"), name="rows")])
def test_axis_types_order_names_and_duplicates_are_preserved(tmp_path, index):
    expected = pd.DataFrame(np.arange(8).reshape(4, 2), index=index, columns=pd.Index(["x", "x"], name="measure"))
    registry, metadata = save_table(tmp_path, expected)
    pd.testing.assert_frame_equal(expected, registry.get(metadata.uri)[1], check_exact=True)


@pytest.mark.parametrize("frame", [pd.DataFrame(index=pd.Index(["001", "NA"], name="id")),
                                  pd.DataFrame({"f": pd.Series(dtype="float32"), "s": pd.Series(dtype="string")}),
                                  pd.DataFrame({"x": [1]}).set_flags(allows_duplicate_labels=False)])
def test_empty_dimensions_and_label_policy_are_preserved(tmp_path, frame):
    registry, metadata = save_table(tmp_path, frame)
    pd.testing.assert_frame_equal(frame, registry.get(metadata.uri)[1], check_exact=True)


@pytest.mark.parametrize("frame", [pd.DataFrame({"x": [[1], [2]]}), pd.DataFrame({"x": [object()]}),
                                  pd.DataFrame({"x": pd.to_datetime(["2026-01-01"])}),
                                  pd.DataFrame({"x": [1]}, index=pd.MultiIndex.from_tuples([("a", 1)]))])
def test_unsupported_semantics_fail_before_immutable_publication(tmp_path, frame):
    storage = ArtifactStorageBackend(str(tmp_path))
    with pytest.raises(ArtifactSerializationError):
        storage.save("table://s/table/v1", frame, ArtifactType.TABLE)
    assert not storage.exists("table://s/table/v1", ArtifactType.TABLE)
    assert not list(tmp_path.rglob("*.tmp"))


def test_legacy_csv_default_decoder_and_version_immutability_are_unchanged(tmp_path):
    storage = ArtifactStorageBackend(str(tmp_path))
    path = tmp_path / "s/table/v1.csv"
    path.parent.mkdir(parents=True)
    path.write_text(",id,mean\n001,NA,12538.393112853733\n002,,0.0\n", encoding="utf-8")
    before = path.read_bytes()
    expected = pd.read_csv(path, index_col=0)
    pd.testing.assert_frame_equal(expected, storage.load("table://s/table/v1", ArtifactType.TABLE), check_exact=True)
    assert expected.iloc[0]["mean"].hex() == "0x1.87d325185a135p+13"
    with pytest.raises(ArtifactAlreadyExistsError):
        storage.save("table://s/table/v1", sample_table(), ArtifactType.TABLE)
    assert path.read_bytes() == before
    assert not (path.parent / "v1.table.npz").exists()


@pytest.mark.parametrize("change", ["dtype", "shape", "decoder", "mask_dtype", "unused_array"])
def test_rewritten_binary_type_contract_fails_closed(tmp_path, change):
    _, metadata = save_table(tmp_path)
    path = Path(metadata.storage_path)
    with np.load(path, allow_pickle=False) as original:
        arrays = {name: np.array(original[name], copy=True) for name in original.files}
    contract = read_contract(path)
    if change == "dtype":
        contract["column_values"][0]["dtype"] = np.dtype("float32").str
    elif change == "shape":
        contract["shape"][0] += 1
    elif change == "decoder":
        contract["decoder"] = "guess_csv_types"
    elif change == "mask_dtype":
        contract["column_values"][6]["mask"]["dtype"] = np.dtype("float64").str
    else:
        arrays["unused"] = np.arange(1)
    arrays[PayloadSerializer.TABLE_CONTRACT_KEY] = np.frombuffer(json.dumps(contract).encode("utf-8"), dtype=np.uint8)
    with path.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    with pytest.raises(ArtifactSerializationError):
        PayloadSerializer.deserialize(path, ArtifactType.TABLE)


def published_table(root):
    registry = ArtifactRegistry(str(root))
    with TaskArtifactTransaction(registry) as transaction:
        transaction.register("table://s/table/v1", sample_table(), ArtifactType.TABLE, "s", "task", "compute")
        transaction.register("json://s/sibling/v1", {"rows": 4}, ArtifactType.JSON, "s", "task", "summarize")
        result = TaskResult(task_id="task", capability="compute", method_used="fixture", status=TaskStatus.SUCCESS,
                            output_artifacts=["table://s/table/v1", "json://s/sibling/v1"])
        registry.commit_task(transaction, "sig", result)
    contract = TaskContract(task_id="task", capability="compute", expected_outputs=result.output_artifacts)
    registry.begin_audit("sig", contract, result)
    registry.record_audit("sig", contract, result, ValidationReport(auditor_name="fixture", target_task_id="task", overall_passed=True))
    return registry, contract, result


def test_typed_table_cross_process_transaction_audit_and_portable_roundtrip(tmp_path):
    registry, contract, result = published_table(tmp_path / "original")
    assert_exact_table(sample_table(), registry.get_audited("table://s/table/v1", "sig", contract)[1])
    bundle = export_artifact_bundle(registry.storage.base_dir, tmp_path / "table.zip")
    destination = import_artifact_bundle(bundle, tmp_path / "relocated")
    reopened = ArtifactRegistry(str(destination / "artifacts"))
    assert_exact_table(sample_table(), reopened.get_audited("table://s/table/v1", "sig", contract)[1])
    assert reopened.get_audit_record("sig") == registry.get_audit_record("sig")
    assert reopened.get_task_commit("sig") == registry.get_task_commit("sig")
    child = tmp_path / "child.py"
    child.write_text("from eacbp.artifact.registry import ArtifactRegistry\n"
                     "from eacbp.schemas.task import TaskContract\n"
                     "import numpy as np, sys\n"
                     "r=ArtifactRegistry(sys.argv[1])\n"
                     "c=TaskContract.model_validate_json(sys.argv[2])\n"
                     "m,t=r.get_audited('table://s/table/v1','sig',c)\n"
                     "assert t.loc['001','mean'].hex()=='0x1.87d325185a136p+13'\n"
                     "assert t['large_int'].iloc[0]==9007199254740993\n"
                     "assert str(t.float32.dtype)=='float32'\n"
                     "assert t.text.tolist()==['001','NA','','NULL']\n"
                     "assert t.nullable_text.isna().tolist()==[False,True,False,False]\n", encoding="utf-8")
    result = subprocess.run([sys.executable, "-B", str(child), str(reopened.storage.base_dir), contract.model_dump_json()],
                            cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    sibling = reopened.get_metadata("json://s/sibling/v1")
    Path(sibling.storage_path).write_text("tampered", encoding="utf-8")
    with pytest.raises(ArtifactAuditAccessError):
        reopened.get_audited("table://s/table/v1", "sig", contract)


def test_unsupported_transaction_sibling_never_commits_partial_task(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    with pytest.raises(ArtifactSerializationError):
        with TaskArtifactTransaction(registry) as transaction:
            transaction.register("table://s/a/v1", sample_table(), ArtifactType.TABLE, "s", "task", "compute")
            transaction.register("table://s/b/v1", pd.DataFrame({"unsupported": [[1]]}), ArtifactType.TABLE, "s", "task", "compute")
    restored = ArtifactRegistry(str(tmp_path))
    assert restored.list_artifacts() == []
    assert restored.get_task_commit("sig") is None
