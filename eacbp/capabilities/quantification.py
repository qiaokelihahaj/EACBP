'''
FASTQ Quantification Capability for EACBP.

The capability has two deliberately separate execution modes:

* ``mode=real`` (the default) requires a complete FASTQ manifest, explicit
  sample metadata, an implemented quantifier, references, and successful
  output for every sample.  It never falls back to generated data.
* ``mode=demo`` is the only mode in which the synthetic SCData generators may
  be used.  Synthetic provenance is copied into both the artifact metadata
  and the SCData ``uns`` payload so downstream consumers cannot mistake it for
  a real quantification.

The input manifest accepts either a pair of paths or lists of paths for each
sample.  Lists are paired by lane and every pair is passed to the quantifier;
the old ``[0]`` behaviour intentionally no longer exists.
'''

from __future__ import annotations

import os
import hashlib
import gzip
import re
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import (
    ExecutionFailureType,
    TaskContract,
    TaskResult,
    TaskStatus,
)


class FASTQManifestError(ValueError):
    """Raised when the registered FASTQ manifest is incomplete or ambiguous."""


def _sha256_file(path: str, chunk_size: int = 1024 * 1024) -> str:
    """Hash a file without loading a sequencing file into memory."""

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


def _safe_component(value: Any, label: str) -> str:
    """Validate a value before it is interpolated into an output path."""

    text = str(value or "")
    if text in {".", ".."} or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", text):
        raise FASTQManifestError(
            f"{label} must be a single safe path component, got {text!r}"
        )
    return text


def _reference_paths(*values: Any) -> List[str]:
    """Expand explicit reference files/directories into stable file paths."""

    paths: List[str] = []
    for value in values:
        if not value:
            continue
        path = Path(value)
        if path.is_file():
            paths.append(str(path))
        elif path.is_dir():
            paths.extend(str(item) for item in sorted(path.rglob("*")) if item.is_file())
    return sorted(set(paths))


def _read_lines(path: Path) -> List[str]:
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\r\n") for line in handle if line.strip()]


def _load_starsolo_matrix(matrix_dir: Path) -> Tuple[SCData, str]:
    """Load a STARsolo Matrix Market directory into the project SCData type."""

    matrix_path = next(
        (path for path in (matrix_dir / "matrix.mtx", matrix_dir / "matrix.mtx.gz") if path.is_file()),
        None,
    )
    barcode_path = next(
        (path for path in (matrix_dir / "barcodes.tsv", matrix_dir / "barcodes.tsv.gz") if path.is_file()),
        None,
    )
    feature_path = next(
        (
            path
            for path in (
                matrix_dir / "features.tsv",
                matrix_dir / "features.tsv.gz",
                matrix_dir / "genes.tsv",
                matrix_dir / "genes.tsv.gz",
            )
            if path.is_file()
        ),
        None,
    )
    if not matrix_path or not barcode_path or not feature_path:
        raise FileNotFoundError(
            f"STARsolo matrix directory is incomplete: {matrix_dir}"
        )

    from scipy.io import mmread

    matrix = mmread(str(matrix_path))
    barcodes = _read_lines(barcode_path)
    features = _read_lines(feature_path)
    if matrix.shape == (len(features), len(barcodes)):
        expression = matrix.T.tocsr() if hasattr(matrix.T, "tocsr") else np.asarray(matrix.T)
    else:
        raise ValueError(
            f"STARsolo matrix shape {matrix.shape} does not match {len(features)} features and {len(barcodes)} barcodes"
        )
    values = expression.data if hasattr(expression, "tocsr") else expression
    if min(expression.shape) == 0 or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("STARsolo counts must be nonempty, finite and nonnegative")
    gene_names = [line.split("\t")[1] if "\t" in line else line.split()[0] for line in features]
    gene_ids = [line.split("\t")[0] for line in features]
    obs = pd.DataFrame(index=pd.Index(barcodes, dtype=str))
    if len(set(gene_ids)) != len(gene_ids) or len(set(barcodes)) != len(barcodes):
        raise ValueError("STARsolo feature IDs and barcodes must be unique")
    var = pd.DataFrame({"gene_name": gene_names}, index=pd.Index(gene_ids, dtype=str))
    data = SCData(expression, obs=obs, var=var)
    data.obs["cell_id"] = data.obs.index.astype(str)
    data.obs["n_counts"] = np.asarray(expression.sum(axis=1)).ravel()
    data.obs["n_genes"] = np.asarray((expression > 0).sum(axis=1)).ravel()
    data.var["n_cells"] = np.asarray((expression > 0).sum(axis=0)).ravel()
    return data, str(matrix_path)


def _as_path_list(value: Any, field: str, sample: str) -> List[str]:
    """Normalize one R1/R2 field without silently dropping lanes."""

    if isinstance(value, (str, os.PathLike)):
        values = [os.fspath(value)]
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray, str)):
        values = [os.fspath(item) for item in value]
    else:
        raise FASTQManifestError(
            f"sample '{sample}' must provide {field} as a path or a list of paths"
        )

    values = [str(path) for path in values if str(path).strip()]
    if not values:
        raise FASTQManifestError(f"sample '{sample}' has no {field} FASTQ files")
    return values


def _lane_key(path: str) -> Tuple[Any, ...]:
    """Return a stable lane key used to pair R1/R2 files."""

    name = Path(path).name
    lane = re.search(r"(?:^|[_-])L(\d{3})(?:[_-]|\.)", name, re.IGNORECASE)
    if lane:
        return (0, int(lane.group(1)))
    # Common read suffixes are removed so R1 and R2 have the same key.
    stem = re.sub(
        r"(?:_R?[12](?:_\d+)?)(?:\.fastq|\.fq)(?:\.gz)?$",
        "",
        name,
        flags=re.IGNORECASE,
    )
    stem = re.sub(r"(?:[_-])R?[12](?:[_-]|$)", "_", stem, flags=re.IGNORECASE)
    return (1, stem.lower())


def _pair_fastqs(r1: Sequence[str], r2: Sequence[str], sample: str) -> List[Tuple[str, str]]:
    """Pair all lane files and reject ambiguous or incomplete manifests."""

    if len(r1) != len(r2):
        raise FASTQManifestError(
            f"sample '{sample}' has {len(r1)} R1 files but {len(r2)} R2 files"
        )

    r1_sorted = sorted((str(p) for p in r1), key=lambda p: (_lane_key(p), p))
    r2_sorted = sorted((str(p) for p in r2), key=lambda p: (_lane_key(p), p))
    keys1 = [_lane_key(p) for p in r1_sorted]
    keys2 = [_lane_key(p) for p in r2_sorted]
    if keys1 != keys2:
        raise FASTQManifestError(
            f"sample '{sample}' R1/R2 lane keys do not match: {keys1!r} vs {keys2!r}"
        )
    return list(zip(r1_sorted, r2_sorted))


def _sample_metadata(
    payload: Mapping[str, Any], sample: str, spec: Mapping[str, Any]
) -> Dict[str, Any]:
    """Read metadata explicitly supplied by the manifest.

    Both ``samples[name].metadata`` and a top-level ``sample_metadata[name]``
    mapping are accepted.  No donor, condition, or batch is inferred from a
    sample name.
    """

    top_level = payload.get("sample_metadata", {})
    external = top_level.get(sample, {}) if isinstance(top_level, Mapping) else {}
    nested = spec.get("metadata", {}) if isinstance(spec, Mapping) else {}
    if nested is None:
        nested = {}
    if external is None:
        external = {}
    if not isinstance(nested, Mapping) or not isinstance(external, Mapping):
        raise FASTQManifestError(f"sample '{sample}' metadata must be a mapping")

    # Top-level fields are accepted for compact manifests, but metadata wins
    # when both representations are present.
    merged = {
        key: value
        for key, value in spec.items()
        if key in {"donor", "condition", "batch", "biological_unit", "timepoint"}
    }
    merged.update(dict(external))
    merged.update(dict(nested))
    missing = [
        key
        for key in ("donor", "condition", "batch")
        if not str(merged.get(key, "")).strip()
    ]
    if missing:
        raise FASTQManifestError(
            f"sample '{sample}' is missing explicit metadata: {', '.join(missing)}"
        )
    merged["sample_id"] = sample
    return merged


def _normalize_samples(payload: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    """Validate sample entries and return all lane pairs plus metadata."""

    raw_samples = payload.get("samples")
    if raw_samples is None:
        raw_samples = payload.get("fastq_paths")
    if not isinstance(raw_samples, Mapping) or not raw_samples:
        raise FASTQManifestError("FASTQ manifest must contain a non-empty 'samples' mapping")

    normalized: Dict[str, Dict[str, Any]] = {}
    for sample_name, raw_spec in raw_samples.items():
        sample = str(sample_name).strip()
        if not sample:
            raise FASTQManifestError("FASTQ manifest contains an empty sample identifier")
        if not isinstance(raw_spec, Mapping):
            raise FASTQManifestError(f"sample '{sample}' specification must be a mapping")

        if raw_spec.get("lanes") is not None:
            lanes = raw_spec.get("lanes")
            if (
                not isinstance(lanes, Sequence)
                or isinstance(lanes, (str, bytes, bytearray))
                or not lanes
            ):
                raise FASTQManifestError(f"sample '{sample}' lanes must be a non-empty list")
            r1: List[str] = []
            r2: List[str] = []
            for index, lane in enumerate(lanes, start=1):
                if not isinstance(lane, Mapping):
                    raise FASTQManifestError(f"sample '{sample}' lane {index} must be a mapping")
                r1_value = lane.get("R1", lane.get("r1"))
                r2_value = lane.get("R2", lane.get("r2"))
                if not r1_value or not r2_value:
                    raise FASTQManifestError(f"sample '{sample}' lane {index} requires R1 and R2")
                r1.append(str(r1_value))
                r2.append(str(r2_value))
            pairs = _pair_fastqs(r1, r2, sample)
        else:
            r1 = _as_path_list(raw_spec.get("R1", raw_spec.get("r1")), "R1", sample)
            r2 = _as_path_list(raw_spec.get("R2", raw_spec.get("r2")), "R2", sample)
            pairs = _pair_fastqs(r1, r2, sample)

        normalized[sample] = {
            "lane_pairs": pairs,
            "metadata": _sample_metadata(payload, sample, raw_spec),
        }
    return normalized


class FASTQQuantificationCapability(BaseCapability):
    """Quantify paired-end FASTQ reads with explicit real/demo semantics."""

    def __init__(
        self,
        capability_name: str = "quantification",
        implementation_id: str = "kb_python_v1",
    ):
        super().__init__(
            capability_name=capability_name,
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_modalities=["scRNA", "FASTQ"],
            accepts_types=[ArtifactType.FASTQ, ArtifactType.JSON],
            suitable_for=[
                "quantification",
                "read_alignment",
                "cell_barcode_demultiplexing",
                "umi_deduplication",
            ],
            output_types=[ArtifactType.ANNDATA],
        )

    def _failure(
        self,
        contract: TaskContract,
        in_uri: Optional[str],
        message: str,
        *,
        status: TaskStatus = TaskStatus.EXECUTION_FAILURE,
        error_type: ExecutionFailureType = ExecutionFailureType.CODE_ERROR,
        metrics: Optional[Dict[str, Any]] = None,
    ) -> TaskResult:
        failure_metrics = dict(metrics or {})
        # A failed request produced no usable data, even when the caller had
        # requested demo mode.  Do not let a failed demo look like a successful
        # synthetic artifact downstream.
        failure_metrics["is_simulated"] = False
        failure_metrics["data_origin"] = "unavailable"
        failure_metrics.setdefault("failure_reason", message)
        return TaskResult(
            task_id=contract.task_id,
            status=status,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri] if in_uri else list(contract.input_artifacts),
            output_artifacts=[],
            metrics=failure_metrics,
            error_type=error_type,
            error_message=message,
        )

    @staticmethod
    def _set_provenance(
        sc_data: SCData,
        *,
        is_simulated: bool,
        data_origin: str,
        extra: Optional[Dict[str, Any]] = None,
    ) -> None:
        provenance = {
            "is_simulated": bool(is_simulated),
            "data_origin": data_origin,
        }
        if extra:
            provenance.update(extra)
        sc_data.uns["quantification_provenance"] = provenance
        sc_data.uns["is_simulated"] = bool(is_simulated)
        sc_data.uns["data_origin"] = data_origin

    def _demo_result(
        self,
        contract: TaskContract,
        registry: ArtifactRegistry,
        in_uri: str,
        study_id: str,
        target_gene: str,
        metrics: Dict[str, Any],
    ) -> TaskResult:
        n_cells = int(contract.parameters.get("n_cells") or 1500)
        n_genes = int(contract.parameters.get("n_genes") or 400)
        random_seed = int(contract.parameters.get("random_seed", 42))
        if "kat8" in target_gene.lower() or "kat8" in study_id.lower() or "cko" in study_id.lower():
            sc_data = SCData.create_synthetic_kat8_study(
                n_cells=n_cells,
                n_genes=n_genes,
                n_cko_mice=4,
                n_con_mice=4,
                random_seed=random_seed,
            )
        else:
            sc_data = SCData.create_synthetic_ad_study(
                n_cells=n_cells,
                n_genes=n_genes,
                n_ad_mice=6,
                n_ctrl_mice=6,
                random_seed=random_seed,
            )

        metrics.update(
            {
                "quant_engine": "native_sc_quantifier_demo",
                "simulated_read_count": int(n_cells * 25000),
                "is_simulated": True,
                "data_origin": "synthetic_demo",
                "samples_quantified": [],
                "samples_failed": {},
                "n_samples_quantified": 0,
            }
        )
        self._set_provenance(sc_data, is_simulated=True, data_origin="synthetic_demo")
        return self._register_success(
            contract,
            registry,
            in_uri,
            study_id,
            sc_data,
            metrics,
            method_used="synthetic_counts_demo_v1",
            executed_operations=["generate_synthetic_counts"],
        )

    def _register_success(
        self,
        contract: TaskContract,
        registry: ArtifactRegistry,
        in_uri: str,
        study_id: str,
        sc_data: SCData,
        metrics: Dict[str, Any],
        *,
        method_used: Optional[str] = None,
        executed_operations: Optional[List[str]] = None,
    ) -> TaskResult:
        metrics.setdefault("n_cells_quantified", sc_data.n_obs)
        metrics.setdefault("n_genes_detected", sc_data.n_vars)
        count_mean = float(sc_data.obs["n_counts"].mean()) if "n_counts" in sc_data.obs.columns else None
        if metrics.get("is_simulated"):
            metrics.setdefault("mean_simulated_counts_per_cell", count_mean)
        else:
            metrics.setdefault("mean_reads_per_cell", count_mean)
        out_uri = contract.expected_outputs[0] if contract.expected_outputs else f"adata://{study_id}/quantified/v1"
        registry.register(
            uri_str=out_uri,
            payload=sc_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=study_id,
            created_by_task=contract.task_id,
            operation="quantify_fastq_reads",
            parent_uris=[in_uri],
            parameters={
                "chemistry": contract.parameters.get("chemistry", "10xv3"),
                "species": contract.parameters.get("species", "unknown"),
                "implementation": self.implementation_id,
                "execution_implementation": method_used or self.implementation_id,
                "mode": contract.parameters.get("mode", "real"),
            },
            summary_metrics=metrics,
        )
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=method_used or self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=executed_operations or [
                "kb_count_alignment",
                "umi_deduplication",
                "gene_annotation_mapping",
            ],
            metrics=metrics,
        )

    def _execute_starsolo_real(
        self,
        contract: TaskContract,
        registry: ArtifactRegistry,
        in_uri: str,
        study_id: str,
        samples: Dict[str, Dict[str, Any]],
        chemistry: str,
        metrics: Dict[str, Any],
    ) -> TaskResult:
        """Run STARsolo for every sample and import its Matrix Market output.

        STARsolo expects cDNA before the barcode/UMI read.  A conventional
        10x manifest calls those files R2 and R1 respectively, so the command
        intentionally emits ``R2 lanes, R1 lanes``.  Each lane is comma
        separated as documented by STARsolo.
        """

        path_env = os.environ.get("PATH", "")
        star_bin = contract.parameters.get("star_bin")
        if not star_bin:
            star_bin = shutil.which("STAR", path=path_env) or shutil.which("star", path=path_env)
        if not star_bin:
            return self._failure(
                contract,
                in_uri,
                "STARsolo requires the STAR executable (set star_bin or add STAR to PATH)",
                error_type=ExecutionFailureType.DEPENDENCY_ERROR,
                metrics=metrics,
            )
        genome_dir_value = contract.parameters.get("genome_dir") or contract.parameters.get("star_genome_dir")
        whitelist_value = contract.parameters.get("whitelist_path") or contract.parameters.get("solo_cb_whitelist")
        genome_dir = Path(genome_dir_value) if genome_dir_value else None
        whitelist_path = Path(whitelist_value) if whitelist_value else None
        missing_refs = []
        if not genome_dir or not genome_dir.is_dir():
            missing_refs.append("genome_dir")
        if not whitelist_path or not whitelist_path.is_file():
            missing_refs.append("whitelist_path")
        if missing_refs:
            return self._failure(
                contract,
                in_uri,
                "STARsolo requires existing references: " + ", ".join(missing_refs),
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )

        gtf_value = contract.parameters.get("gtf_path") or contract.parameters.get("sjdb_gtf")
        gtf_path = Path(gtf_value) if gtf_value else None
        if gtf_value and not gtf_path.is_file():
            return self._failure(
                contract,
                in_uri,
                f"STARsolo gtf_path does not exist: {gtf_value}",
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )

        try:
            input_paths = sorted(
                {
                    path
                    for spec in samples.values()
                    for pair in spec["lane_pairs"]
                    for path in pair
                }
            )
            reference_paths = _reference_paths(str(genome_dir), str(whitelist_path), str(gtf_path) if gtf_path else None)
            if not reference_paths:
                raise OSError("STARsolo references contain no files")
            input_hashes = {path: _sha256_file(path) for path in input_paths}
            reference_hashes = {path: _sha256_file(path) for path in reference_paths}
        except OSError as exc:
            return self._failure(
                contract,
                in_uri,
                f"failed to hash STARsolo inputs or references: {exc}",
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )

        metrics.update(
            {
                "input_file_hashes": input_hashes,
                "reference_hashes": reference_hashes,
                "commands": [],
            }
        )
        try:
            run_id = _safe_component(contract.parameters.get("run_id") or uuid.uuid4().hex, "run_id")
            _safe_component(study_id, "study_id")
        except FASTQManifestError as exc:
            return self._failure(contract, in_uri, str(exc), metrics=metrics)
        base_work_dir = Path(contract.parameters.get("work_dir") or f"outputs/starsolo_{study_id}")
        out_dir = base_work_dir / f"run_{run_id}"
        try:
            out_dir.mkdir(parents=True, exist_ok=False)
        except OSError as exc:
            return self._failure(
                contract,
                in_uri,
                f"cannot create isolated STARsolo run directory '{out_dir}': {exc}",
                error_type=ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT,
                metrics=metrics,
            )

        threads = int(contract.parameters.get("threads", 8))
        umi_len = "12" if chemistry == "10xv3" else "10"
        sample_adatas: List[Any] = []
        successful_samples: List[str] = []
        failed_samples: Dict[str, str] = {}

        try:
            import anndata as ad
        except Exception as exc:
            return self._failure(
                contract,
                in_uri,
                f"anndata is required to load STARsolo output: {exc}",
                error_type=ExecutionFailureType.DEPENDENCY_ERROR,
                metrics=metrics,
            )

        for sample_name, sample_spec in samples.items():
            try:
                _safe_component(sample_name, "sample name")
            except FASTQManifestError as exc:
                failed_samples[sample_name] = str(exc)
                continue
            sample_out = out_dir / sample_name
            sample_out.mkdir(parents=True, exist_ok=False)
            lane_pairs = sample_spec["lane_pairs"]
            all_reads = [path for pair in lane_pairs for path in pair]
            compressed = [path.lower().endswith(".gz") for path in all_reads]
            if any(compressed) and not all(compressed):
                failed_samples[sample_name] = "STARsolo requires all lanes to use the same compression mode"
                continue
            command = [
                str(star_bin),
                "--runThreadN",
                str(threads),
                "--genomeDir",
                str(genome_dir),
                "--readFilesIn",
                ",".join(r2 for _r1, r2 in lane_pairs),
                ",".join(r1 for r1, _r2 in lane_pairs),
                "--soloType",
                "CB_UMI_Simple",
                "--soloCBwhitelist",
                str(whitelist_path),
                "--soloCBstart",
                "1",
                "--soloCBlen",
                "16",
                "--soloUMIstart",
                "17",
                "--soloUMIlen",
                umi_len,
                "--soloFeatures",
                "Gene",
                "--outFileNamePrefix",
                str(sample_out) + os.sep,
            ]
            if all(compressed):
                command.extend(["--readFilesCommand", "zcat"])
            if gtf_path:
                command.extend(["--sjdbGTFfile", str(gtf_path)])
            metrics["commands"].append(command)
            try:
                subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    check=True,
                    env={**os.environ, "PATH": path_env},
                )
            except (OSError, subprocess.CalledProcessError) as exc:
                stderr = getattr(exc, "stderr", "") or str(exc)
                failed_samples[sample_name] = str(stderr)[:1000]
                continue

            matrix_root = sample_out / "Solo.out" / "Gene"
            matrix_candidates = [matrix_root / ("raw" if contract.parameters.get("use_raw_counts", False) else "filtered")]
            loaded = None
            load_error = None
            for matrix_dir in matrix_candidates:
                try:
                    loaded = _load_starsolo_matrix(matrix_dir)
                    break
                except (FileNotFoundError, ValueError, OSError) as exc:
                    load_error = exc
            if loaded is None:
                failed_samples[sample_name] = f"STARsolo matrix output could not be loaded: {load_error}"
                continue
            sample_data, matrix_path = loaded
            metadata = sample_spec["metadata"]
            for key, value in metadata.items():
                if key != "sample_id":
                    sample_data.obs[key] = value
            sample_data.obs["sample"] = sample_name
            sample_data.obs["sample_id"] = sample_name
            sample_data.obs["source_fastq_lanes"] = ",".join(
                f"{Path(r1).name}|{Path(r2).name}" for r1, r2 in lane_pairs
            )
            sample_data.obs.index = [f"{sample_name}_{barcode}" for barcode in sample_data.obs.index]
            sample_data.obs["cell_id"] = sample_data.obs.index
            sample_adatas.append(sample_data.to_anndata())
            successful_samples.append(sample_name)

        metrics["samples_quantified"] = successful_samples
        metrics["samples_failed"] = failed_samples
        metrics["n_samples_quantified"] = len(successful_samples)
        current_hashes = {}
        try:
            current_hashes.update({path: _sha256_file(path) for path in input_paths})
            current_hashes.update({path: _sha256_file(path) for path in reference_paths})
        except OSError as exc:
            return self._failure(
                contract,
                in_uri,
                f"failed to re-hash STARsolo inputs or references: {exc}",
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )
        expected_hashes = {**input_hashes, **reference_hashes}
        if current_hashes != expected_hashes:
            return self._failure(
                contract,
                in_uri,
                "FASTQ or STARsolo reference contents changed during quantification",
                error_type=ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT,
                metrics=metrics,
            )
        if failed_samples or len(successful_samples) != len(samples):
            return self._failure(
                contract,
                in_uri,
                "STARsolo quantification failed for one or more samples",
                error_type=ExecutionFailureType.CODE_ERROR,
                metrics=metrics,
            )
        if not sample_adatas:
            return self._failure(
                contract,
                in_uri,
                "STARsolo produced no sample matrices",
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )
        combined = sample_adatas[0] if len(sample_adatas) == 1 else ad.concat(sample_adatas, join="outer", fill_value=0.0)
        sc_data = SCData.from_anndata(combined)
        self._set_provenance(
            sc_data,
            is_simulated=False,
            data_origin="real_fastq",
            extra={
                "samples_quantified": successful_samples,
                "lane_pairs_per_sample": {name: len(spec["lane_pairs"]) for name, spec in samples.items()},
                "input_file_hashes": input_hashes,
                "reference_hashes": reference_hashes,
            },
        )
        metrics.update({"quant_engine": "STARsolo", "is_simulated": False, "data_origin": "real_fastq"})
        return self._register_success(
            contract,
            registry,
            in_uri,
            study_id,
            sc_data,
            metrics,
            method_used=self.implementation_id,
            executed_operations=["starsolo_alignment", "umi_deduplication", "gene_annotation_mapping"],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0] if contract.input_artifacts else None
        mode = str(contract.parameters.get("mode", "real")).strip().lower()
        if mode not in {"real", "demo"}:
            return self._failure(
                contract,
                in_uri,
                f"unsupported quantification mode '{mode}'; expected 'real' or 'demo'",
            )
        if not in_uri:
            return self._failure(contract, in_uri, "FASTQ quantification requires an input artifact URI")

        try:
            _meta, payload = registry.get(in_uri)
        except Exception as exc:
            return self._failure(contract, in_uri, f"FASTQ input artifact is unavailable: {exc}")

        if not isinstance(payload, Mapping):
            return self._failure(contract, in_uri, "FASTQ input artifact payload must be a manifest mapping")

        parsed_uri = ArtifactURI.parse(in_uri)
        study_id = contract.parameters.get("study_id") or parsed_uri.study_id
        try:
            _safe_component(study_id, "study_id")
        except FASTQManifestError as exc:
            return self._failure(contract, in_uri, str(exc))
        target_gene = str(contract.parameters.get("target_gene", "Kat8"))
        metrics: Dict[str, Any] = {
            "mode": mode,
            "is_simulated": mode == "demo",
            "data_origin": "synthetic_demo" if mode == "demo" else "real_fastq",
        }

        # STARsolo demo is intentionally not represented as a STARsolo run;
        # callers must use real mode when they request this executable method.
        if self.implementation_id.lower().startswith("starsolo") and mode == "demo":
            return self._failure(
                contract,
                in_uri,
                "STARsolo demo execution is unavailable; use mode=real with STAR and references",
                status=TaskStatus.METHOD_FAILURE,
                error_type=ExecutionFailureType.CODE_ERROR,
                metrics=metrics,
            )

        if mode == "demo":
            # Demo still validates a registered manifest when one is present,
            # but does not require files or external executables.  This makes
            # synthetic tests explicit while keeping their setup lightweight.
            try:
                normalized = _normalize_samples(payload)
                metrics["samples_requested"] = list(normalized)
                metrics["lane_pairs_requested"] = sum(len(v["lane_pairs"]) for v in normalized.values())
            except FASTQManifestError:
                # A demo with no reads is allowed by design; it remains clearly
                # synthetic and cannot be mistaken for a real run.
                metrics["samples_requested"] = []
                metrics["lane_pairs_requested"] = 0
            return self._demo_result(contract, registry, in_uri, study_id, target_gene, metrics)

        # From this point onward execution is real-only.  Any problem returns a
        # structured failure and no output artifact is registered.
        try:
            samples = _normalize_samples(payload)
        except FASTQManifestError as exc:
            return self._failure(
                contract,
                in_uri,
                str(exc),
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )

        metrics["samples_requested"] = list(samples)
        metrics["lane_pairs_requested"] = sum(len(v["lane_pairs"]) for v in samples.values())
        missing_files = [
            path
            for sample in samples.values()
            for pair in sample["lane_pairs"]
            for path in pair
            if not Path(path).is_file()
        ]
        if missing_files:
            return self._failure(
                contract,
                in_uri,
                "FASTQ input files are missing: " + ", ".join(missing_files[:10]),
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )

        if self.implementation_id.lower().startswith("starsolo"):
            chemistry_value = contract.parameters.get("chemistry")
            if chemistry_value is None:
                chemistry_value = payload.get("chemistry")
            chemistry = str(chemistry_value or "").strip()
            if chemistry not in {"10xv2", "10xv3"}:
                return self._failure(
                    contract,
                    in_uri,
                    "real FASTQ quantification requires explicit supported chemistry ('10xv2' or '10xv3')",
                    error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                    metrics=metrics,
                )
            return self._execute_starsolo_real(
                contract,
                registry,
                in_uri,
                study_id,
                samples,
                chemistry,
                metrics,
            )

        if self.implementation_id != "kb_python_v1":
            return self._failure(
                contract,
                in_uri,
                f"real FASTQ quantification implementation '{self.implementation_id}' is unavailable",
                status=TaskStatus.METHOD_FAILURE,
                error_type=ExecutionFailureType.DEPENDENCY_ERROR,
                metrics=metrics,
            )

        py_bin_dir = str(Path(sys.executable).parent)
        path_env = os.environ.get("PATH", "")
        if py_bin_dir not in path_env.split(os.pathsep):
            path_env = py_bin_dir + os.pathsep + path_env
        kb_bin = shutil.which("kb", path=path_env)
        if not kb_bin:
            candidate = Path(sys.executable).parent / "kb"
            kb_bin = str(candidate) if candidate.is_file() else None
        if not kb_bin:
            return self._failure(
                contract,
                in_uri,
                "kb-python executable 'kb' is unavailable for real FASTQ quantification",
                error_type=ExecutionFailureType.DEPENDENCY_ERROR,
                metrics=metrics,
            )

        # References are intentionally explicit.  The capability does not
        # consult historical or machine-specific default paths.
        index_path = contract.parameters.get("index_path")
        t2g_path = contract.parameters.get("t2g_path")
        missing_refs = [
            name
            for name, value in (("index_path", index_path), ("t2g_path", t2g_path))
            if not value or not Path(value).is_file()
        ]
        if missing_refs:
            return self._failure(
                contract,
                in_uri,
                "real FASTQ quantification requires existing references: " + ", ".join(missing_refs),
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )

        run_id = str(contract.parameters.get("run_id") or uuid.uuid4().hex)
        base_work_dir = Path(contract.parameters.get("work_dir") or f"outputs/kb_quant_{study_id}")
        # Every invocation gets a private subdirectory.  Existing files are
        # never treated as a successful cache.
        out_dir = base_work_dir / f"run_{run_id}"
        try:
            out_dir.mkdir(parents=True, exist_ok=False)
        except OSError as exc:
            return self._failure(
                contract,
                in_uri,
                f"cannot create isolated quantification run directory '{out_dir}': {exc}",
                error_type=ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT,
                metrics=metrics,
            )
        chemistry_value = contract.parameters.get("chemistry")
        if chemistry_value is None:
            chemistry_value = payload.get("chemistry")
        chemistry = str(chemistry_value or "").strip()
        if chemistry not in {"10xv2", "10xv3"}:
            return self._failure(
                contract,
                in_uri,
                "real FASTQ quantification requires explicit supported chemistry ('10xv2' or '10xv3')",
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )
        threads = int(contract.parameters.get("threads", 8))
        num_reads = contract.parameters.get("num_reads")
        sample_adatas: List[Any] = []
        successful_samples: List[str] = []
        failed_samples: Dict[str, str] = {}

        try:
            import anndata as ad
        except Exception as exc:
            return self._failure(
                contract,
                in_uri,
                f"anndata is required to load real kb-python output: {exc}",
                error_type=ExecutionFailureType.DEPENDENCY_ERROR,
                metrics=metrics,
            )

        input_paths = sorted(
            {
                path
                for spec in samples.values()
                for pair in spec["lane_pairs"]
                for path in pair
            }
        )
        reference_paths = [str(index_path), str(t2g_path)]
        try:
            input_hashes = {path: _sha256_file(path) for path in input_paths}
            reference_hashes = {path: _sha256_file(path) for path in reference_paths}
        except OSError as exc:
            return self._failure(
                contract,
                in_uri,
                f"failed to hash FASTQ inputs or kb-python references: {exc}",
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )
        metrics["input_file_hashes"] = input_hashes
        metrics["reference_hashes"] = reference_hashes
        metrics["commands"] = []
        for sample_name, sample_spec in samples.items():
            sample_out = out_dir / sample_name
            sample_out.mkdir(parents=True, exist_ok=False)
            expected_h5ad = sample_out / "counts_unfiltered" / "adata.h5ad"
            lane_pairs = sample_spec["lane_pairs"]
            command = [
                str(kb_bin),
                "count",
                "-i",
                str(index_path),
                "-g",
                str(t2g_path),
                "-x",
                chemistry,
                "-o",
                str(sample_out),
                "--h5ad",
                "--gene-names",
                "-t",
                str(threads),
                "-m",
                "8G",
            ]
            if num_reads:
                command.extend(["-N", str(num_reads)])
            # Preserve every lane pair in input order.  No historical output
            # is loaded before this command succeeds.
            command.extend([path for pair in lane_pairs for path in pair])
            metrics["commands"].append(command)
            try:
                subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    check=True,
                    env={**os.environ, "PATH": path_env},
                )
            except (OSError, subprocess.CalledProcessError) as exc:
                stderr = getattr(exc, "stderr", "") or str(exc)
                failed_samples[sample_name] = str(stderr)[:1000]
                continue

            if not expected_h5ad.is_file():
                failed_samples[sample_name] = "kb-python completed without creating counts_unfiltered/adata.h5ad"
                continue
            try:
                sample_adata = ad.read_h5ad(str(expected_h5ad))
                metadata = sample_spec["metadata"]
                for key, value in metadata.items():
                    if key != "sample_id":
                        sample_adata.obs[key] = value
                sample_adata.obs["sample"] = sample_name
                sample_adata.obs["sample_id"] = sample_name
                sample_adata.obs["source_fastq_lanes"] = ",".join(
                    f"{Path(r1).name}|{Path(r2).name}" for r1, r2 in lane_pairs
                )
                sample_adata.obs_names = [f"{sample_name}_{barcode}" for barcode in sample_adata.obs_names]
                sample_adatas.append(sample_adata)
                successful_samples.append(sample_name)
            except Exception as exc:
                failed_samples[sample_name] = f"failed to load kb-python output: {exc}"

        metrics["samples_quantified"] = successful_samples
        metrics["samples_failed"] = failed_samples
        metrics["n_samples_quantified"] = len(successful_samples)
        try:
            current_input_hashes = {path: _sha256_file(path) for path in input_paths}
            current_reference_hashes = {path: _sha256_file(path) for path in reference_paths}
        except OSError as exc:
            return self._failure(
                contract,
                in_uri,
                f"failed to re-hash FASTQ inputs or kb-python references: {exc}",
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )
        if current_input_hashes != input_hashes or current_reference_hashes != reference_hashes:
            return self._failure(
                contract,
                in_uri,
                "FASTQ or kb-python reference contents changed during quantification",
                error_type=ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT,
                metrics=metrics,
            )
        if failed_samples or len(successful_samples) != len(samples):
            return self._failure(
                contract,
                in_uri,
                "real FASTQ quantification failed for one or more samples",
                error_type=ExecutionFailureType.CODE_ERROR,
                metrics=metrics,
            )
        if not sample_adatas:
            return self._failure(
                contract,
                in_uri,
                "real FASTQ quantification produced no sample matrices",
                error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                metrics=metrics,
            )

        if len(sample_adatas) == 1:
            combined_adata = sample_adatas[0]
        else:
            combined_adata = ad.concat(sample_adatas, join="outer", fill_value=0.0)
        sc_data = SCData.from_anndata(combined_adata)
        self._set_provenance(
            sc_data,
            is_simulated=False,
            data_origin="real_fastq",
            extra={
                "samples_quantified": successful_samples,
                "lane_pairs_per_sample": {
                    name: len(spec["lane_pairs"]) for name, spec in samples.items()
                },
                "input_file_hashes": input_hashes,
                "reference_hashes": reference_hashes,
            },
        )
        metrics.update(
            {
                "quant_engine": "kb-python (kallisto | bustools)",
                "is_simulated": False,
                "data_origin": "real_fastq",
            }
        )
        return self._register_success(contract, registry, in_uri, study_id, sc_data, metrics)
