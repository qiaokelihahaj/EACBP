from __future__ import annotations

import hashlib
import io
import json
from dataclasses import replace
from pathlib import Path, PureWindowsPath

import pytest

from scripts.fetch_public_validation_data import (
    DATASETS,
    DownloadError,
    _download_one,
    _receipt_entry,
    _write_receipt,
    validate_root,
)


def _spec(relative_path: str, content: bytes):
    return replace(
        DATASETS[0],
        key="test",
        relative_path=relative_path,
        url="https://example.invalid/data.h5ad",
        sha256=hashlib.sha256(content).hexdigest(),
        md5=None,
    )


class _Response(io.BytesIO):
    pass


class _Opener:
    def __init__(self, payload: bytes) -> None:
        self.payload = payload
        self.calls = 0

    def open(self, request, timeout):
        self.calls += 1
        return _Response(self.payload)


def test_root_requires_absolute_external_directory_and_rejects_windows_c_drive(tmp_path: Path) -> None:
    repository = tmp_path / "checkout"
    repository.mkdir()

    with pytest.raises(ValueError, match="absolute"):
        validate_root("relative/data", repository_root=repository, windows_platform=False)
    with pytest.raises(ValueError, match="outside"):
        validate_root(str(repository / "data"), repository_root=repository, windows_platform=False)
    with pytest.raises(ValueError, match="C:"):
        validate_root(
            r"C:\EACBP-public-validation",
            repository_root=r"C:\Users\user\EACBP",
            windows_platform=True,
        )
    with pytest.raises(ValueError, match="C:"):
        validate_root(
            r"\\?\C:\EACBP-public-validation",
            repository_root=r"C:\Users\user\EACBP",
            windows_platform=True,
        )

    accepted = validate_root(
        r"F:\EACBP-public-validation",
        repository_root=r"C:\Users\user\EACBP",
        windows_platform=True,
    )
    assert PureWindowsPath(accepted).drive == "F:"


def test_dataset_catalog_records_verified_license_provenance_without_inventing_derivative_terms() -> None:
    raw, processed, kang = DATASETS

    assert raw.license_name == "CC BY 4.0 (upstream 10x PBMC3k dataset)"
    assert processed.license_url == raw.license_url
    assert "does not publish a separate license" in processed.license_note
    assert processed.provenance_url is not None
    assert kang.license_name == "CC BY 4.0"
    assert kang.source_version == "2"
    assert kang.source_file_id == "34464122"
    assert kang.md5 == "adb2246232e8493031c576982c0c02a3"
    assert raw.relative_path == "datasets/pbmc3k_raw.h5ad"
    assert processed.relative_path == "references/pbmc3k_processed.h5ad"
    assert kang.relative_path == "datasets/kang_2018.h5ad"


def test_existing_valid_file_is_verified_and_skipped(tmp_path: Path) -> None:
    content = b"verified h5ad placeholder"
    spec = _spec("datasets/sample.h5ad", content)
    target = tmp_path / spec.relative_path
    target.parent.mkdir()
    target.write_bytes(content)
    opener = _Opener(b"must not be used")

    result = _download_one(spec, tmp_path, opener, timeout=1)

    assert result["downloaded"] is False
    assert result["size_bytes"] == len(content)
    assert opener.calls == 0
    assert target.read_bytes() == content


def test_download_verifies_then_atomically_replaces_bad_file(tmp_path: Path) -> None:
    expected = b"complete verified payload"
    spec = _spec("datasets/sample.h5ad", expected)
    target = tmp_path / spec.relative_path
    target.parent.mkdir()
    target.write_bytes(b"old contents")

    result = _download_one(spec, tmp_path, _Opener(expected), timeout=1)

    assert result["downloaded"] is True
    assert result["sha256"] == spec.sha256
    assert target.read_bytes() == expected
    assert not (target.parent / (target.name + ".part")).exists()


def test_failed_checksum_keeps_previous_target_and_cleans_partial(tmp_path: Path) -> None:
    previous = b"pre-existing target"
    spec = _spec("datasets/sample.h5ad", b"expected but not received")
    target = tmp_path / spec.relative_path
    target.parent.mkdir()
    target.write_bytes(previous)

    with pytest.raises(DownloadError, match="checksum mismatch"):
        _download_one(spec, tmp_path, _Opener(b"truncated response"), timeout=1)

    assert target.read_bytes() == previous
    assert not (target.parent / (target.name + ".part")).exists()


def test_receipt_records_source_license_checksums_and_size(tmp_path: Path) -> None:
    content = b"receipt payload"
    spec = _spec("datasets/sample.h5ad", content)
    result = {
        "relative_path": spec.relative_path,
        "size_bytes": len(content),
        "sha256": spec.sha256,
        "md5": hashlib.md5(content).hexdigest(),
        "downloaded": True,
    }

    receipt_path = _write_receipt(tmp_path, [_receipt_entry(spec, result)])

    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    dataset = receipt["datasets"][0]
    assert dataset["download_url"] == spec.url
    assert dataset["source"]["url"] == spec.source_url
    assert dataset["citation"] == spec.citation
    assert dataset["license"]["name"] == spec.license_name
    assert dataset["checksums"]["sha256"] == spec.sha256
    assert dataset["size_bytes"] == len(content)
