#!/usr/bin/env python3
"""Fetch checksum-pinned public datasets for validation into external storage."""

from __future__ import annotations

import argparse
import hashlib
import json
import ntpath
import os
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath
from typing import Any
from urllib.parse import urlsplit


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = "F:/EACBP-public-validation" if os.name == "nt" else "/mnt/f/EACBP-public-validation"
RECEIPT_NAME = "receipt.json"
BUFFER_SIZE = 1024 * 1024

TENX_PBMC3K_URL = "https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0"
TENX_PBMC3K_LICENSE_URL = "https://creativecommons.org/licenses/by/4.0/"
SCANPY_RAW_URL = "https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.datasets.pbmc3k.html"
SCANPY_PROCESSED_URL = "https://scanpy.readthedocs.io/en/1.12.x/generated/scanpy.datasets.pbmc3k_processed.html"
KANG_FIGSHARE_URL = "https://figshare.com/articles/dataset/Kang_HM_Subramaniam_M_Targ_S_Nguyen_M_et_al_2017/19397624"
KANG_LICENSE_URL = "https://creativecommons.org/licenses/by/4.0/"


@dataclass(frozen=True)
class DatasetSpec:
    key: str
    relative_path: str
    url: str
    sha256: str
    source_title: str
    source_url: str
    provenance_url: str | None
    citation: str
    license_name: str
    license_url: str | None
    license_note: str
    md5: str | None = None
    source_version: str | None = None
    source_file_id: str | None = None


DATASETS = (
    DatasetSpec(
        key="pbmc3k_raw",
        relative_path="datasets/pbmc3k_raw.h5ad",
        url="https://exampledata.scverse.org/scanpy/pbmc3k_raw.h5ad",
        sha256="89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1",
        source_title="10x Genomics 3k PBMCs from a Healthy Donor; Scanpy raw serialization",
        source_url=TENX_PBMC3K_URL,
        provenance_url=SCANPY_RAW_URL,
        citation=(
            "Zheng, G.X.Y. et al. Massively parallel digital transcriptional profiling of "
            "single cells. Nature Communications 8, 14049 (2017). doi:10.1038/ncomms14049. "
            "Dataset: 10x Genomics, 3k PBMCs from a Healthy Donor."
        ),
        license_name="CC BY 4.0 (upstream 10x PBMC3k dataset)",
        license_url=TENX_PBMC3K_LICENSE_URL,
        license_note=(
            "The official 10x PBMC3k page states CC BY 4.0. Scanpy describes this H5AD as "
            "a serialization of that source data; it does not publish a separate file-specific license."
        ),
    ),
    DatasetSpec(
        key="pbmc3k_processed",
        relative_path="references/pbmc3k_processed.h5ad",
        url="https://exampledata.scverse.org/scanpy/pbmc3k.h5ad",
        sha256="0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38",
        source_title="10x Genomics 3k PBMCs from a Healthy Donor; Scanpy processed reference",
        source_url=TENX_PBMC3K_URL,
        provenance_url=SCANPY_PROCESSED_URL,
        citation=(
            "Wolf, F.A., Angerer, P. & Theis, F.J. SCANPY: large-scale single-cell gene "
            "expression data analysis. Genome Biology 19, 15 (2018). "
            "doi:10.1186/s13059-017-1382-0. Dataset: 10x Genomics, 3k PBMCs from a Healthy Donor."
        ),
        license_name="CC BY 4.0 (upstream 10x PBMC3k dataset)",
        license_url=TENX_PBMC3K_LICENSE_URL,
        license_note=(
            "The official 10x PBMC3k page states CC BY 4.0 for the source dataset. Scanpy says "
            "this H5AD is tutorial-processed but does not publish a separate license for this derivative."
        ),
    ),
    DatasetSpec(
        key="kang_2018",
        relative_path="datasets/kang_2018.h5ad",
        url="https://exampledata.scverse.org/pertpy/kang_2018.h5ad",
        sha256="e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830",
        md5="adb2246232e8493031c576982c0c02a3",
        source_title="Kang HM, Subramaniam M, Targ S, Nguyen M et al. 2017 (Figshare version 2)",
        source_url=KANG_FIGSHARE_URL,
        provenance_url=None,
        source_version="2",
        source_file_id="34464122",
        citation=(
            "Kang, H.M. et al. Multiplexed droplet single-cell RNA-sequencing using natural "
            "genetic variation. Nature Biotechnology (2018). doi:10.1038/nbt.4042. "
            "Figshare dataset DOI:10.6084/m9.figshare.19397624, version 2."
        ),
        license_name="CC BY 4.0",
        license_url=KANG_LICENSE_URL,
        license_note="The Figshare article lists CC BY 4.0 for version 2.",
    ),
)


class DownloadError(RuntimeError):
    """Raised when a downloaded file does not match its pinned checksums."""


def _is_relative_to(path: Any, parent: Any) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _is_windows_c_drive(path: PureWindowsPath) -> bool:
    drive = path.drive.upper()
    return drive == "C:" or drive in {"\\\\?\\C:", "\\\\.\\C:"}


def validate_root(
    root_value: str | os.PathLike[str],
    *,
    repository_root: str | os.PathLike[str] = REPOSITORY_ROOT,
    windows_platform: bool | None = None,
) -> str:
    """Validate an absolute data root outside the checkout and return its normalized path."""
    use_windows_paths = os.name == "nt" if windows_platform is None else windows_platform
    value = os.fspath(root_value)
    repository = os.fspath(repository_root)

    if use_windows_paths:
        root = PureWindowsPath(ntpath.normpath(value))
        if not root.is_absolute():
            raise ValueError("--root must be an absolute Windows path")
        if _is_windows_c_drive(root):
            raise ValueError("--root on Windows must not use the C: drive")
        repo = PureWindowsPath(ntpath.normpath(repository))
        if _is_relative_to(root, repo):
            raise ValueError("--root must be outside the repository checkout")
        if os.name == "nt":
            resolved = Path(str(root)).resolve(strict=False)
            resolved_windows = PureWindowsPath(str(resolved))
            if _is_windows_c_drive(resolved_windows):
                raise ValueError("--root resolves to the C: drive, which is not allowed")
            if _is_relative_to(resolved_windows, repo):
                raise ValueError("--root resolves inside the repository checkout")
            return str(resolved_windows)
        return str(root)

    root_path = Path(value).expanduser()
    if not root_path.is_absolute():
        raise ValueError("--root must be an absolute path")
    root_path = root_path.resolve(strict=False)
    repo_path = Path(repository).expanduser().resolve(strict=False)
    if _is_relative_to(root_path, repo_path):
        raise ValueError("--root must be outside the repository checkout")
    return str(root_path)


def _build_opener(proxy: str | None) -> urllib.request.OpenerDirector:
    if proxy is None:
        return urllib.request.build_opener()
    parsed = urlsplit(proxy)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        raise ValueError("--proxy must be an absolute http:// or https:// URL")
    return urllib.request.build_opener(
        urllib.request.ProxyHandler({"http": proxy, "https": proxy})
    )


def _hash_file(path: Path) -> tuple[str, str, int]:
    sha256 = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False)
    size_bytes = 0
    with path.open("rb") as stream:
        while chunk := stream.read(BUFFER_SIZE):
            sha256.update(chunk)
            md5.update(chunk)
            size_bytes += len(chunk)
    return sha256.hexdigest(), md5.hexdigest(), size_bytes


def _checksums_match(spec: DatasetSpec, sha256: str, md5: str) -> bool:
    return sha256 == spec.sha256 and (spec.md5 is None or md5 == spec.md5)


def _download_one(
    spec: DatasetSpec,
    root: Path,
    opener: urllib.request.OpenerDirector,
    *,
    timeout: float,
) -> dict[str, int | str | bool]:
    relative_path = Path(spec.relative_path)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise DownloadError(f"dataset path must stay under the selected root: {spec.relative_path}")
    target = root / relative_path
    if target.is_symlink():
        raise DownloadError(f"refusing to replace symlink target: {target}")
    if target.exists():
        if not target.is_file():
            raise DownloadError(f"target exists and is not a regular file: {target}")
        old_sha256, old_md5, old_size = _hash_file(target)
        if _checksums_match(spec, old_sha256, old_md5):
            return {
                "relative_path": relative_path.as_posix(),
                "size_bytes": old_size,
                "sha256": old_sha256,
                "md5": old_md5,
                "downloaded": False,
            }
        print(
            f"Existing {target.name} failed checksum validation; fetching a verified replacement.",
            file=sys.stderr,
        )

    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".part")
    if part.is_symlink():
        raise DownloadError(f"refusing to write through symlink partial file: {part}")
    try:
        request = urllib.request.Request(
            spec.url,
            headers={"User-Agent": "EACBP-public-validation-data/1.0"},
        )
        with opener.open(request, timeout=timeout) as response, part.open("wb") as output:
            sha256 = hashlib.sha256()
            md5 = hashlib.md5(usedforsecurity=False)
            size_bytes = 0
            while chunk := response.read(BUFFER_SIZE):
                output.write(chunk)
                sha256.update(chunk)
                md5.update(chunk)
                size_bytes += len(chunk)
            output.flush()
            os.fsync(output.fileno())

        observed_sha256 = sha256.hexdigest()
        observed_md5 = md5.hexdigest()
        if not _checksums_match(spec, observed_sha256, observed_md5):
            expected = f"SHA-256 {spec.sha256}"
            observed = f"SHA-256 {observed_sha256}"
            if spec.md5:
                expected += f", MD5 {spec.md5}"
                observed += f", MD5 {observed_md5}"
            raise DownloadError(
                f"checksum mismatch for {spec.relative_path}: expected {expected}; got {observed}"
            )

        os.replace(part, target)
        return {
            "relative_path": relative_path.as_posix(),
            "size_bytes": size_bytes,
            "sha256": observed_sha256,
            "md5": observed_md5,
            "downloaded": True,
        }
    except BaseException:
        try:
            part.unlink(missing_ok=True)
        except OSError:
            pass
        raise


def _receipt_entry(spec: DatasetSpec, result: dict[str, int | str | bool]) -> dict[str, Any]:
    source: dict[str, str] = {"title": spec.source_title, "url": spec.source_url}
    if spec.source_version is not None:
        source["version"] = spec.source_version
    if spec.source_file_id is not None:
        source["file_id"] = spec.source_file_id
    if spec.provenance_url is not None:
        source["processing_or_provenance_url"] = spec.provenance_url
    checksums: dict[str, str] = {"sha256": spec.sha256}
    if spec.md5 is not None:
        checksums["md5"] = spec.md5
    return {
        "key": spec.key,
        "relative_path": result["relative_path"],
        "download_url": spec.url,
        "source": source,
        "citation": spec.citation,
        "license": {
            "name": spec.license_name,
            "url": spec.license_url,
            "note": spec.license_note,
        },
        "checksums": checksums,
        "verified_sha256": result["sha256"],
        "verified_md5": result["md5"],
        "size_bytes": result["size_bytes"],
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def _write_receipt(root: Path, entries: list[dict[str, Any]]) -> Path:
    receipt = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "datasets": entries,
    }
    target = root / RECEIPT_NAME
    if target.is_symlink():
        raise DownloadError(f"refusing to replace symlink receipt: {target}")
    part = target.with_name(target.name + ".part")
    if part.is_symlink():
        raise DownloadError(f"refusing to write through symlink receipt partial: {part}")
    try:
        with part.open("w", encoding="utf-8", newline="\n") as stream:
            json.dump(receipt, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(part, target)
    except BaseException:
        try:
            part.unlink(missing_ok=True)
        except OSError:
            pass
        raise
    return target


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download the pinned PBMC3k raw/processed references and Kang 2018 H5AD into "
            "an absolute directory outside this checkout. Existing valid files are verified and skipped."
        ),
        epilog=(
            "Example (PowerShell): python scripts/fetch_public_validation_data.py "
            "--root F:/EACBP-public-validation --proxy http://127.0.0.1:7890"
        ),
    )
    parser.add_argument("--root", default=DEFAULT_ROOT, help=f"external absolute data directory (default: {DEFAULT_ROOT})")
    parser.add_argument("--proxy", help="optional HTTP(S) proxy URL, for example http://127.0.0.1:7890")
    parser.add_argument("--timeout", type=float, default=120.0, help="per-request timeout in seconds (default: 120)")
    args = parser.parse_args(argv)
    if args.timeout <= 0:
        parser.error("--timeout must be greater than zero")
    return args


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        root_text = validate_root(args.root)
        root = Path(root_text).resolve(strict=False)
        opener = _build_opener(args.proxy)
        entries = []
        for spec in DATASETS:
            result = _download_one(spec, root, opener, timeout=args.timeout)
            action = "downloaded" if result["downloaded"] else "verified existing"
            print(
                f"{action}: {result['relative_path']} ({result['size_bytes']} bytes, "
                f"SHA-256 {result['sha256']})"
            )
            entries.append(_receipt_entry(spec, result))
        receipt_path = _write_receipt(root, entries)
        print(f"receipt: {receipt_path}")
        return 0
    except (OSError, urllib.error.URLError, ValueError, DownloadError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
