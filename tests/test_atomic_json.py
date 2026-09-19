"""Failure-path coverage for the shared atomic JSON writer."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from eacbp._atomic_json import atomic_write_json


@pytest.mark.parametrize("failure", ["dump", "fsync", "replace"])
def test_atomic_json_failure_preserves_previous_file_and_cleans_stage(tmp_path, failure):
    target = tmp_path / "document.json"
    target.write_text('{"old": 1}\n', encoding="utf-8")
    before = target.read_bytes()

    if failure == "dump":
        patch_target, error = "eacbp._atomic_json.json.dump", RuntimeError("dump failed")
    elif failure == "fsync":
        patch_target, error = "eacbp._atomic_json.os.fsync", OSError("fsync failed")
    else:
        patch_target, error = "eacbp._atomic_json.os.replace", OSError("replace failed")

    with patch(patch_target, side_effect=error), pytest.raises(type(error), match=str(error)):
        atomic_write_json(target, {"new": 2}, dump_kwargs={"allow_nan": False})

    assert target.read_bytes() == before
    assert list(tmp_path.glob(f".{target.name}.*.tmp")) == []
