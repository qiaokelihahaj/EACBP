"""Optional real single-cell QC, annotation, and background-removal methods.

The classes in this module intentionally do not implement a lightweight
fallback.  Scrublet and CellTypist are called through their local Python APIs,
and CellBender is called through an explicitly supplied executable.  This is
important for provenance: a missing optional dependency or external program
must be visible as a failed task instead of looking like a successful
analysis performed by a different algorithm.

The three public capabilities are:

``ScrubletDoubletCapability``
    ``doublet_detection / scanpy_scrublet_v1``.  Scores raw counts separately
    for each library batch, marks doublets by default, and filters only when
    ``filter_doublets=True`` is explicit.

``CellTypistAnnotationCapability``
    ``cell_annotation / celltypist_local_v1``.  Loads a local model path only;
    it never asks CellTypist to download a named model.  Existing cell type
    annotations are retained unless ``use_as_cell_type=True`` is explicit.

``CellBenderBackgroundRemovalCapability``
    ``background_removal / cellbender_cli_v1``.  Runs a user-supplied
    CellBender executable against a user-supplied unfiltered droplet input and
    stores the resulting matrix in ``layers['corrected_counts']``.

The module is deliberately independent of the capability factory.  The
factory/planner can register these classes when the optional dependencies and
external-resource parameters are available.
"""

from __future__ import annotations

import copy
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _coerce_scdata(payload: Any) -> SCData:
    """Return an independent :class:`SCData` copy for a registry payload."""

    if isinstance(payload, SCData):
        return payload.copy()
    if hasattr(payload, "X") and hasattr(payload, "obs") and hasattr(payload, "var"):
        return SCData.from_anndata(payload)
    if isinstance(payload, Mapping):
        return SCData.from_dict(dict(payload))
    raise TypeError(
        "Advanced QC capabilities require an AnnData/SCData artifact; "
        f"received {type(payload).__name__}"
    )


def _load_scdata(contract: TaskContract, registry: ArtifactRegistry) -> Tuple[str, SCData]:
    if not contract.input_artifacts:
        raise ValueError("Advanced QC capabilities require one input AnnData artifact")
    in_uri = ArtifactURI.parse(contract.input_artifacts[0]).to_string()
    _, payload = registry.get(in_uri)
    return in_uri, _coerce_scdata(payload)


def _study_id(in_uri: str) -> str:
    return ArtifactURI.parse(in_uri).study_id


def _primary_output_uri(
    contract: TaskContract,
    in_uri: str,
    default_path: str,
) -> str:
    """Resolve the first output while keeping expected_outputs authoritative."""

    if contract.expected_outputs:
        uri = ArtifactURI.parse(contract.expected_outputs[0]).to_string()
    else:
        uri = f"adata://{_study_id(in_uri)}/{default_path}"
    parsed = ArtifactURI.parse(uri)
    if parsed.scheme != "adata":
        raise ValueError(
            "The first expected output of an advanced QC capability must be "
            f"an AnnData URI, received {uri!r}"
        )
    if parsed.study_id != _study_id(in_uri):
        raise ValueError("Output study_id must match the input artifact study_id")
    return uri


def _optional_output_uri(contract: TaskContract, index: int) -> Optional[str]:
    if len(contract.expected_outputs) <= index:
        return None
    if len(contract.expected_outputs) > 2:
        raise ValueError("Advanced QC capabilities support at most two expected outputs")
    return ArtifactURI.parse(contract.expected_outputs[index]).to_string()


def _output_type(uri: str) -> ArtifactType:
    scheme = ArtifactURI.parse(uri).scheme
    if scheme == "table":
        return ArtifactType.TABLE
    if scheme == "json":
        return ArtifactType.JSON
    raise ValueError("Optional advanced QC output must use table:// or json://")


def _register_outputs(
    *,
    contract: TaskContract,
    registry: ArtifactRegistry,
    in_uri: str,
    payload: SCData,
    capability: BaseCapability,
    primary_default_path: str,
    operation: str,
    metrics: Dict[str, Any],
    software_versions: Dict[str, str],
    optional_payload: Any = None,
    optional_operation: Optional[str] = None,
) -> TaskResult:
    """Register an AnnData output and an optional table/JSON audit output."""

    primary_uri = _primary_output_uri(contract, in_uri, primary_default_path)
    # Validate the complete output contract before publishing the primary
    # artifact.  This prevents a malformed second URI from leaving a partial
    # successful-looking result in a direct capability invocation.
    optional_uri = _optional_output_uri(contract, 1)
    if optional_uri is not None:
        if ArtifactURI.parse(optional_uri).study_id != _study_id(in_uri):
            raise ValueError("Optional output study_id must match the input artifact study_id")
        optional_type = _output_type(optional_uri)
        if optional_payload is None:
            raise ValueError("A second expected output was supplied but no optional payload was produced")
    sid = _study_id(in_uri)
    registry.register(
        uri_str=primary_uri,
        payload=payload.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id=sid,
        created_by_task=contract.task_id,
        operation=operation,
        parent_uris=[in_uri],
        parameters=contract.parameters,
        software_versions=software_versions,
        random_seed=int(contract.parameters.get("random_seed", 42)),
        summary_metrics=metrics,
    )

    outputs = [primary_uri]
    if optional_uri is not None:
        registry.register(
            uri_str=optional_uri,
            payload=optional_payload,
            artifact_type=optional_type,
            study_id=sid,
            created_by_task=contract.task_id,
            operation=optional_operation or f"{operation}_audit",
            parent_uris=[primary_uri],
            parameters=contract.parameters,
            software_versions=software_versions,
            random_seed=int(contract.parameters.get("random_seed", 42)),
            summary_metrics=metrics,
        )
        outputs.append(optional_uri)

    return TaskResult(
        task_id=contract.task_id,
        status=TaskStatus.SUCCESS,
        capability=capability.capability_name,
        method_used=capability.implementation_id,
        input_artifacts=[in_uri],
        output_artifacts=outputs,
        executed_operations=list(getattr(capability, "contract_operations", [])),
        metrics=metrics,
    )


def _package_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def _matrix_values(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.data)
    return np.asarray(matrix)


def _copy_matrix(matrix: Any) -> Any:
    if sparse.issparse(matrix):
        return matrix.copy()
    return np.array(matrix, copy=True)


def _ensure_counts_matrix(matrix: Any, label: str) -> None:
    """Validate that a matrix can be used as raw UMI/count input."""

    values = _matrix_values(matrix)
    if values.size == 0:
        raise ValueError(f"{label} is empty")
    if not np.isfinite(values).all():
        raise ValueError(f"{label} contains non-finite values")
    if (values < 0).any():
        raise ValueError(f"{label} contains negative values")
    # Raw counts can be stored as float32 by AnnData, but fractional values are
    # not counts and should never be silently rounded for Scrublet/CellBender.
    if not np.allclose(values, np.rint(values), rtol=0.0, atol=1e-6):
        raise ValueError(f"{label} must contain integer-valued raw counts")


def _counts_for_scrublet(data: SCData, params: Mapping[str, Any]) -> Tuple[Any, str]:
    layer_name = params.get("counts_layer", "counts")
    if layer_name is not None and str(layer_name) in data.layers:
        name = str(layer_name)
        counts = data.layers[name]
        source = f"layer:{name}"
    else:
        counts = data.X
        source = "X"
    _ensure_counts_matrix(counts, f"Scrublet counts ({source})")
    if getattr(counts, "shape", None) != data.shape:
        raise ValueError(f"Scrublet counts ({source}) shape does not match the input AnnData")
    return counts, source


def _normalize_hash(value: Any) -> str:
    value = str(value).strip().lower()
    return value.removeprefix("sha256:")


def _sha256_path(path: Path) -> str:
    """Hash a file or a directory deterministically using SHA-256.

    Directory inputs are common for unfiltered 10x matrices.  Their relative
    file names are included so that renaming a matrix file cannot pass an
    integrity check accidentally.
    """

    digest = hashlib.sha256()
    if path.is_file():
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    if path.is_dir():
        files = sorted(item for item in path.rglob("*") if item.is_file())
        if not files:
            raise ValueError(f"External input directory is empty: {path}")
        for item in files:
            relative = item.relative_to(path).as_posix().encode("utf-8")
            digest.update(len(relative).to_bytes(8, "big"))
            digest.update(relative)
            with item.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
    else:
        raise FileNotFoundError(f"External resource does not exist: {path}")
    return digest.hexdigest()


def _verify_external_hashes(params: Mapping[str, Any], observed: Mapping[str, str]) -> None:
    expected = params.get("external_resource_sha256")
    if expected is None:
        return
    if not isinstance(expected, Mapping):
        raise ValueError("external_resource_sha256 must be a mapping of parameter name to SHA-256")
    for key, value in expected.items():
        key = str(key)
        actual = observed.get(key)
        if actual is None:
            # Permit planners to key the map by the resolved path, while the
            # normal contract uses parameter names such as model_path.
            actual = observed.get(str(Path(key).resolve()))
        if actual is None:
            raise ValueError(f"No observed external resource hash is available for {key!r}")
        if _normalize_hash(actual) != _normalize_hash(value):
            raise ValueError(
                f"External resource hash mismatch for {key!r}: "
                f"expected {_normalize_hash(value)}, observed {_normalize_hash(actual)}"
            )


def _as_json_value(value: Any) -> Any:
    """Convert model/CLI metadata to a conservative JSON-compatible value."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (np.generic,)):
        return value.item()
    if isinstance(value, Mapping):
        return {str(k): _as_json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_as_json_value(v) for v in value]
    return str(value)


def _cellbender_sidecar_paths(output_path: Path, params: Mapping[str, Any]) -> Dict[str, Path]:
    """Resolve CellBender's optional report, metrics and log sidecars.

    CellBender writes ``<output stem>_report.html`` and
    ``<output stem>_metrics.csv``.  The log name has varied between releases
    and invocation wrappers, so the adapter accepts an explicit path and
    otherwise checks the common ``.log``/``_log`` spellings.  Missing sidecars
    are recorded by the caller rather than silently treated as quality proof.
    """

    candidates = {
        "report": [
            output_path.with_name(f"{output_path.stem}_report.html"),
            output_path.with_name(f"{output_path.stem}.html"),
        ],
        "metrics": [
            output_path.with_name(f"{output_path.stem}_metrics.csv"),
            output_path.with_name(f"{output_path.stem}.csv"),
        ],
        "log": [
            output_path.with_name(f"{output_path.stem}.log"),
            output_path.with_name(f"{output_path.stem}_log"),
            output_path.with_name(f"{output_path.stem}_log.txt"),
        ],
    }
    resolved: Dict[str, Path] = {}
    for kind, values in candidates.items():
        explicit = params.get(f"{kind}_path")
        if explicit:
            resolved[kind] = Path(str(explicit)).expanduser().resolve()
            continue
        existing = next((candidate for candidate in values if candidate.exists()), None)
        resolved[kind] = (existing or values[0]).resolve()
    return resolved


def _cellbender_checkpoint_resources(extra_args: Sequence[str], run_cwd: Path) -> List[Tuple[str, Path]]:
    """Resolve local ``--checkpoint`` arguments for external hash checking."""

    resources: List[Tuple[str, Path]] = []
    index = 0
    while index < len(extra_args):
        argument = extra_args[index]
        checkpoint_value = None
        if argument == "--checkpoint":
            if index + 1 >= len(extra_args):
                raise ValueError("CellBender extra_args/cli_args has --checkpoint without a path")
            checkpoint_value = extra_args[index + 1]
            index += 2
        elif argument.startswith("--checkpoint="):
            checkpoint_value = argument.split("=", 1)[1]
            index += 1
        else:
            index += 1
        if checkpoint_value is None:
            continue
        if not checkpoint_value.strip():
            raise ValueError("CellBender extra_args/cli_args has an empty --checkpoint path")
        checkpoint = Path(checkpoint_value).expanduser()
        if not checkpoint.is_absolute():
            checkpoint = run_cwd / checkpoint
        checkpoint = checkpoint.resolve(strict=True)
        if not checkpoint.is_file():
            raise ValueError(f"CellBender checkpoint must be a file: {checkpoint}")
        key = "checkpoint" if not resources else f"checkpoint[{len(resources)}]"
        resources.append((key, checkpoint))
    return resources


def _read_cellbender_metrics(path: Path) -> Dict[str, Any]:
    """Read CellBender's key/value metrics CSV conservatively.

    CellBender 0.4 writes a two-column CSV without a header.  Supporting a
    header as well keeps the provenance parser compatible with wrapper scripts
    while bounding the serialized payload to scalar values.
    """

    frame = pd.read_csv(path, header=None, dtype=str)
    if frame.empty:
        return {}
    if frame.shape[1] < 2:
        return {"rows": int(len(frame))}
    metrics: Dict[str, Any] = {}
    for key, value in frame.iloc[:, :2].itertuples(index=False, name=None):
        key = str(key).strip()
        if not key or key.casefold() in {"metric", "key", "name"}:
            continue
        raw = str(value).strip()
        try:
            number = float(raw)
            metrics[key] = int(number) if number.is_integer() else number
        except (TypeError, ValueError):
            metrics[key] = raw
        if len(metrics) >= 200:
            break
    return metrics


def _cellbender_quality_warnings(
    *,
    stdout: str,
    stderr: str,
    log_text: str,
    report_text: str,
    sidecar_metrics: Mapping[str, Any],
) -> List[str]:
    """Extract bounded external quality warnings without inventing a pass."""

    # Keep lines rather than entire HTML/log files in the artifact metadata.
    # Console output is intentionally broad because wrappers and the
    # CellBender process report warnings there.  The HTML report needs a
    # separate, narrower matcher after removing non-visible CSS and scripts.
    stream_pattern = re.compile(
        r"warning|warn\b|suboptimal|convergence|jupyter.*not found|could be",
        re.IGNORECASE,
    )
    report_pattern = re.compile(
        r"(?:\bwarning\b\s*:"
        r"|\belbo\b.*(?:\bdeviat(?:es|ed|ing)\b|\bwrong direction\b|\bsuboptimal\b)"
        r"|\b(?:deviat(?:es|ed|ing)\b|\bwrong direction\b|\bsuboptimal\b).*\belbo\b"
        r"|\b(?:output|result)\b.*\bsuboptimal\b"
        r"|\bunusual behavior\b.*\blearning[- ]rate\b)",
        re.IGNORECASE,
    )

    stream_warnings: List[str] = []
    report_warnings: List[str] = []
    metric_warnings: List[str] = []

    def _add_warning(target: List[str], source: str, line: Any) -> None:
        normalized = " ".join(str(line).split())
        # HTML tags such as ``<em>WARNING</em>:`` produce a space before the
        # colon when separate text nodes are joined.  Keep the stored evidence
        # readable and stable across equivalent markup.
        normalized = re.sub(r"\s+([,:;.!?])", r"\1", normalized)
        if not normalized:
            return
        item = f"{source}: {normalized[:500]}"
        if item not in target:
            target.append(item)

    for source, text in (
        ("stdout", stdout),
        ("stderr", stderr),
        ("log", log_text),
    ):
        if not text:
            continue
        for line in str(text).splitlines():
            normalized = " ".join(line.split())
            if normalized and stream_pattern.search(normalized):
                _add_warning(stream_warnings, source, normalized)

    if report_text:
        # CellBender's report is a Jupyter HTML export.  Searching its raw
        # source finds dozens of Jupyter CSS declarations such as
        # ``--jp-warn-color0`` before reaching the actual ELBO commentary.
        # Parse only visible body text and discard style/script content.  This
        # deliberately uses the standard library so report parsing does not
        # add an optional dependency to the adapter.
        from html.parser import HTMLParser

        class _VisibleReportParser(HTMLParser):
            _ignored_tags = {"style", "script", "noscript", "template"}
            _block_tags = {
                "address", "article", "aside", "blockquote", "body", "br",
                "caption", "dd", "div", "dl", "dt", "fieldset", "figcaption",
                "figure", "footer", "form", "h1", "h2", "h3", "h4", "h5",
                "h6", "head", "header", "hr", "li", "main", "nav", "ol",
                "p", "pre", "section", "table", "tbody", "td", "tfoot", "th",
                "thead", "tr", "ul",
            }

            def __init__(self) -> None:
                super().__init__(convert_charrefs=True)
                self.lines: List[str] = []
                self._buffer: List[str] = []
                self._body_seen = False
                self._body_depth = 0
                self._head_depth = 0
                self._ignored_depth = 0

            def _active(self) -> bool:
                return bool(
                    not self._head_depth
                    and not self._ignored_depth
                    and (not self._body_seen or self._body_depth > 0)
                )

            def _flush(self) -> None:
                if not self._buffer:
                    return
                text = " ".join("".join(self._buffer).split())
                self._buffer = []
                text = re.sub(r"\s+([,:;.!?])", r"\1", text)
                if text:
                    self.lines.append(text)

            def handle_starttag(self, tag: str, attrs: Any) -> None:
                tag = tag.casefold()
                if tag == "body":
                    self._flush()
                    self._body_seen = True
                    self._body_depth += 1
                    return
                if tag == "head":
                    self._flush()
                    self._head_depth += 1
                    return
                if tag in self._ignored_tags:
                    self._flush()
                    self._ignored_depth += 1
                    return
                if self._active() and tag in self._block_tags:
                    self._flush()

            def handle_endtag(self, tag: str) -> None:
                tag = tag.casefold()
                if tag in self._ignored_tags and self._ignored_depth:
                    self._ignored_depth -= 1
                    return
                if tag == "head" and self._head_depth:
                    self._head_depth -= 1
                    return
                if tag == "body" and self._body_depth:
                    self._flush()
                    self._body_depth -= 1
                    return
                if self._active() and tag in self._block_tags:
                    self._flush()

            def handle_data(self, data: str) -> None:
                if self._active():
                    self._buffer.append(data)

        parser = _VisibleReportParser()
        try:
            parser.feed(str(report_text))
            parser.close()
            parser._flush()
        except (TypeError, ValueError):
            # A malformed report should not make the external run look like a
            # successful clean-quality run.  Stream and metric evidence still
            # remain available to the caller.
            parser.lines = []
        for line in parser.lines:
            if report_pattern.search(line):
                _add_warning(report_warnings, "report", line)

    # CellBender's metrics CSV can carry a convergence indicator even when the
    # HTML report was not generated.  Preserve it as evidence; do not turn a
    # version-specific numeric threshold into an automatic failure.
    indicator = sidecar_metrics.get("convergence_indicator")
    try:
        if indicator is not None and float(indicator) > 1.0:
            metric_warnings.append(
                "metrics: convergence_indicator="
                f"{float(indicator):g}; inspect the CellBender report before scientific use"
            )
    except (TypeError, ValueError):
        pass

    # The metadata remains bounded, while report/metric evidence is protected
    # from a long stream of ordinary process warnings.  In the usual case this
    # preserves the original stdout -> stderr -> log -> report ordering.  On
    # overflow, stderr/log and report/metric evidence are retained before
    # lower-priority stdout entries so an important ELBO warning cannot be
    # displaced by the cap.
    combined = stream_warnings + report_warnings + metric_warnings
    if len(combined) <= 32:
        return combined
    protected = report_warnings + metric_warnings
    remaining = max(0, 32 - len(protected))
    if len(protected) >= 32:
        return protected[:32]
    stderr_log = [
        item for item in stream_warnings
        if item.startswith("stderr: ") or item.startswith("log: ")
    ]
    stdout_items = [item for item in stream_warnings if item.startswith("stdout: ")]
    selected_stream = stderr_log[:remaining]
    remaining -= len(selected_stream)
    selected_stream.extend(stdout_items[:remaining])
    return selected_stream + protected


def _collect_cellbender_sidecars(
    output_path: Path,
    params: Mapping[str, Any],
    *,
    stdout: str,
    stderr: str,
) -> Dict[str, Any]:
    """Capture sidecar provenance and quality caveats after a CLI run."""

    paths = _cellbender_sidecar_paths(output_path, params)
    sidecar: Dict[str, Any] = {}
    report_text = ""
    log_text = ""
    metrics_values: Dict[str, Any] = {}
    missing: List[str] = []
    for kind, path in paths.items():
        available = path.is_file()
        sidecar[f"{kind}_path"] = str(path)
        sidecar[f"{kind}_available"] = bool(available)
        sidecar[f"{kind}_sha256"] = _sha256_path(path) if available else None
        if not available:
            missing.append(kind)
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except (OSError, UnicodeError):
            text = ""
        if kind == "report":
            report_text = text
        elif kind == "log":
            log_text = text
        elif kind == "metrics":
            try:
                metrics_values = _read_cellbender_metrics(path)
            except (OSError, ValueError, pd.errors.ParserError) as exc:
                sidecar["metrics_parse_error"] = f"{type(exc).__name__}: {exc}"

    warnings = _cellbender_quality_warnings(
        stdout=stdout,
        stderr=stderr,
        log_text=log_text,
        report_text=report_text,
        sidecar_metrics=metrics_values,
    )
    if missing:
        warnings.extend(
            f"sidecar: CellBender {kind} sidecar was not found at {paths[kind]}"
            for kind in missing
        )
    warnings = list(dict.fromkeys(warnings))[:32]
    report_available = bool(sidecar["report_available"])
    if not report_available:
        quality_status = "completed_report_unavailable"
    elif warnings:
        quality_status = "completed_with_quality_warnings"
    else:
        quality_status = "completed_with_report"
    sidecar.update(
        {
            "metrics_values": _as_json_value(metrics_values),
            "quality_warnings": warnings,
            "quality_status": quality_status,
            # This adapter records external evidence only.  A scientific
            # quality gate still needs an independent review of the report and
            # convergence diagnostics.
            "scientific_quality_certified": False,
        }
    )
    return sidecar


def _single_metadata_value(data: SCData, field: str) -> Optional[str]:
    """Read a single species/tissue value from obs or uns when available."""

    values: List[Any] = []
    if field in data.obs.columns:
        values = [v for v in data.obs[field].tolist() if pd.notna(v)]
    elif field in data.uns:
        value = data.uns.get(field)
        if value is not None:
            values = [value]
    if not values:
        return None
    unique = {str(v).strip() for v in values if str(v).strip()}
    return next(iter(unique)) if len(unique) == 1 else None


def _cell_and_gene_names(data: SCData, params: Mapping[str, Any]) -> Tuple[np.ndarray, np.ndarray]:
    cell_key = str(params.get("cell_id_key", "cell_id"))
    gene_key = str(params.get("gene_id_key", "gene_name"))
    if cell_key in data.obs.columns:
        cells = data.obs[cell_key].astype(str).to_numpy()
    else:
        cells = data.obs.index.astype(str).to_numpy()
    if gene_key in data.var.columns:
        genes = data.var[gene_key].astype(str).to_numpy()
    else:
        genes = data.var.index.astype(str).to_numpy()
    if len(set(cells)) != len(cells):
        raise ValueError("Input cell identifiers are not unique; alignment cannot be certified")
    if len(set(genes)) != len(genes):
        raise ValueError("Input gene identifiers are not unique; alignment cannot be certified")
    return cells, genes


def _external_anndata_names(adata: Any, params: Mapping[str, Any]) -> Tuple[np.ndarray, np.ndarray]:
    cell_key = str(params.get("cell_id_key", "cell_id"))
    gene_key = str(params.get("gene_id_key", "gene_name"))
    if cell_key in adata.obs.columns:
        cells = adata.obs[cell_key].astype(str).to_numpy()
    else:
        cells = np.asarray(adata.obs_names.astype(str))
    if gene_key in adata.var.columns:
        genes = adata.var[gene_key].astype(str).to_numpy()
    else:
        genes = np.asarray(adata.var_names.astype(str))
    if len(set(cells)) != len(cells) or len(set(genes)) != len(genes):
        raise ValueError("CellBender output contains duplicate cell or gene identifiers")
    return cells, genes


def _map_external_matrix_to_target(
    matrix: Any,
    external_cells: np.ndarray,
    external_genes: np.ndarray,
    target_cells: np.ndarray,
    target_genes: np.ndarray,
) -> Tuple[Any, Dict[str, Any]]:
    """Safely map an external matrix onto the target artifact's IDs.

    CellBender's full ``output.h5`` contains all input droplet barcodes, while
    the current artifact commonly contains only a filtered population.  A
    dimension/order equality check would therefore reject a valid run.  The
    mapping is deliberately stricter than an intersection: every target ID
    must occur exactly once in the external output, and duplicate IDs are
    rejected before any indexing.  Extra external droplets/features are
    allowed and recorded in the mapping report.
    """

    def build_index(values: np.ndarray, label: str) -> Dict[str, int]:
        strings = [str(value) for value in values]
        if len(strings) != len(set(strings)):
            raise ValueError(f"CellBender external {label} identifiers are duplicated")
        return {value: index for index, value in enumerate(strings)}

    ext_cell_index = build_index(external_cells, "cell")
    ext_gene_index = build_index(external_genes, "gene")
    target_cells = np.asarray([str(value) for value in target_cells], dtype=object)
    target_genes = np.asarray([str(value) for value in target_genes], dtype=object)
    missing_cells = [value for value in target_cells if value not in ext_cell_index]
    missing_genes = [value for value in target_genes if value not in ext_gene_index]
    if missing_cells:
        raise ValueError(
            "CellBender output is missing target cell IDs: "
            f"{missing_cells[:10]}" + (" ..." if len(missing_cells) > 10 else "")
        )
    if missing_genes:
        raise ValueError(
            "CellBender output is missing target gene IDs: "
            f"{missing_genes[:10]}" + (" ..." if len(missing_genes) > 10 else "")
        )

    row_positions = np.asarray([ext_cell_index[value] for value in target_cells], dtype=int)
    col_positions = np.asarray([ext_gene_index[value] for value in target_genes], dtype=int)
    if sparse.issparse(matrix):
        mapped = matrix[row_positions][:, col_positions].copy()
    else:
        mapped = np.asarray(matrix)[np.ix_(row_positions, col_positions)].copy()
    mapping = {
        "target_n_cells": int(len(target_cells)),
        "target_n_genes": int(len(target_genes)),
        "external_n_cells": int(len(external_cells)),
        "external_n_genes": int(len(external_genes)),
        "mapped_n_cells": int(len(row_positions)),
        "mapped_n_genes": int(len(col_positions)),
        "extra_external_cells": int(len(external_cells) - len(target_cells)),
        "extra_external_genes": int(len(external_genes) - len(target_genes)),
        "cell_order_exact": bool(np.array_equal(external_cells, target_cells)),
        "gene_order_exact": bool(np.array_equal(external_genes, target_genes)),
        "cell_ids_unique": True,
        "gene_ids_unique": True,
    }
    return mapped, mapping


def _label_is_unknown(value: Any) -> bool:
    if value is None or pd.isna(value):
        return True
    normalized = str(value).strip().casefold()
    return normalized in {
        "",
        "unknown",
        "unassigned",
        "unclassifiable",
        "undefined",
        "none",
        "nan",
        "na",
        "n/a",
    }


def _clean_subprocess_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


# ---------------------------------------------------------------------------
# Scanpy/Scrublet
# ---------------------------------------------------------------------------


class ScrubletDoubletCapability(BaseCapability):
    """Score doublets with Scanpy's real Scrublet implementation."""

    contract_operations = [
        "validate_raw_counts",
        "score_doublets_by_batch",
        "record_doublet_threshold",
        "mark_predicted_doublets",
    ]

    def resolve_contract_operations(self, parameters):
        operations = list(self.contract_operations)
        if parameters.get("filter_doublets", False):
            operations.append("filter_cells")
        return operations

    def __init__(self):
        super().__init__(
            capability_name="doublet_detection",
            implementation_id="scanpy_scrublet_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA, ArtifactType.TABLE, ArtifactType.JSON],
        )
        self.legacy_aliases = {
            "scrublet": "scanpy_scrublet_v1",
            "scrublet_v1": "scanpy_scrublet_v1",
        }

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri, data = _load_scdata(contract, registry)
        params = contract.parameters
        if data.n_obs < 3 or data.n_vars < 3:
            raise ValueError("Scrublet requires at least three cells and three genes")

        counts, counts_source = _counts_for_scrublet(data, params)
        batch_key = params.get("batch_key", params.get("batch_col", "batch"))
        if batch_key is not None:
            batch_key = str(batch_key)
            if batch_key not in data.obs.columns:
                raise ValueError(
                    f"Scrublet batch_key {batch_key!r} is missing; supply library-batch metadata "
                    "or explicitly set batch_key=None for a single-batch run"
                )
            if data.obs[batch_key].isna().any():
                raise ValueError("Scrublet batch metadata must be complete")

        try:
            scanpy = importlib.import_module("scanpy")
            # Current Scanpy ships its own Scrublet implementation. Requiring
            # the unrelated standalone package adds an unnecessary annoy/C++
            # build dependency on Windows.
            if not callable(getattr(scanpy.pp, "scrublet", None)):
                raise ImportError("Installed Scanpy does not expose pp.scrublet")
        except ImportError as exc:
            raise ImportError(
                "scanpy_scrublet_v1 requires Scanpy with pp.scrublet; no fallback was run"
            ) from exc

        # Work on a separate AnnData object with counts in X.  This allows an
        # upstream normalized X to stay untouched while still enforcing raw
        # counts for the actual Scrublet call.
        working = data.to_anndata()
        working.X = _copy_matrix(counts)
        if "counts" not in working.layers:
            working.layers["counts"] = _copy_matrix(counts)

        requested_n_pc = int(params.get("n_prin_comps", 30))
        max_n_pc = min(working.n_obs - 1, working.n_vars - 1)
        if max_n_pc < 1:
            raise ValueError("Scrublet cannot compute a principal component on this matrix")
        effective_n_pc = min(max(1, requested_n_pc), max_n_pc)

        threshold = params.get("threshold", None)
        if threshold is not None:
            threshold = float(threshold)
            if not np.isfinite(threshold) or threshold < 0:
                raise ValueError("Scrublet threshold must be a finite nonnegative number")

        scrublet_kwargs: Dict[str, Any] = {
            "batch_key": batch_key,
            "sim_doublet_ratio": float(params.get("sim_doublet_ratio", 2.0)),
            "expected_doublet_rate": float(params.get("expected_doublet_rate", 0.05)),
            "stdev_doublet_rate": float(params.get("stdev_doublet_rate", 0.02)),
            "synthetic_doublet_umi_subsampling": float(
                params.get("synthetic_doublet_umi_subsampling", 1.0)
            ),
            "knn_dist_metric": params.get("knn_dist_metric", "euclidean"),
            "normalize_variance": bool(params.get("normalize_variance", True)),
            "log_transform": bool(params.get("log_transform", False)),
            "mean_center": bool(params.get("mean_center", True)),
            "n_prin_comps": effective_n_pc,
            "use_approx_neighbors": params.get("use_approx_neighbors", None),
            "get_doublet_neighbor_parents": bool(params.get("get_doublet_neighbor_parents", False)),
            "threshold": threshold,
            "verbose": bool(params.get("verbose", False)),
            "copy": False,
            "random_state": int(params.get("random_seed", 42)),
        }
        if "n_neighbors" in params and params.get("n_neighbors") is not None:
            scrublet_kwargs["n_neighbors"] = int(params["n_neighbors"])

        try:
            scanpy.pp.scrublet(working, **scrublet_kwargs)
        except (ImportError, ModuleNotFoundError) as exc:
            raise ImportError(
                "scanpy_scrublet_v1 could not import its Scrublet runtime; "
                "install the missing dependency reported by Scanpy and retry"
            ) from exc

        if "doublet_score" not in working.obs.columns:
            raise ValueError("Scanpy Scrublet returned no doublet_score column")
        score = np.asarray(working.obs["doublet_score"].to_numpy(), dtype=float)
        if score.shape != (data.n_obs,) or not np.isfinite(score).all():
            raise ValueError("Scanpy Scrublet returned invalid doublet scores")
        prediction_key = "predicted_doublet"
        if prediction_key not in working.obs.columns:
            prediction_key = "predicted_doublets"
        if prediction_key not in working.obs.columns:
            raise ValueError("Scanpy Scrublet returned no predicted-doublet column")
        predicted = np.asarray(working.obs[prediction_key].to_numpy(), dtype=bool)
        if predicted.shape != (data.n_obs,):
            raise ValueError("Scanpy Scrublet returned predictions with the wrong cell count")

        scrublet_uns = copy.deepcopy(dict(working.uns.get("scrublet", {})))
        threshold_by_batch: Dict[str, Optional[float]] = {}
        if batch_key is not None:
            batches_uns = scrublet_uns.get("batches", {})
            if isinstance(batches_uns, Mapping):
                for batch, batch_record in batches_uns.items():
                    value = batch_record.get("threshold") if isinstance(batch_record, Mapping) else None
                    threshold_by_batch[str(batch)] = float(value) if value is not None else None
        else:
            value = scrublet_uns.get("threshold")
            threshold_by_batch["__all__"] = float(value) if value is not None else None
        thresholds = [value for value in threshold_by_batch.values() if value is not None]

        output = data.copy()
        output.obs["doublet_score"] = score
        output.obs["predicted_doublet"] = predicted
        if "counts" not in output.layers:
            # Make the raw matrix explicit even when the caller supplied it
            # in X instead of a named layer.
            output.layers["counts"] = _copy_matrix(counts)
        output.uns["scrublet"] = scrublet_uns
        output.uns["eacbp_scrublet"] = {
            "implementation_id": self.implementation_id,
            "counts_source": counts_source,
            "counts_layer_preserved": "counts" in output.layers or counts_source == "X",
            "batch_key": batch_key,
            "threshold_requested": threshold,
            "threshold_by_batch": threshold_by_batch,
            "filter_doublets": bool(params.get("filter_doublets", False)),
        }

        filter_doublets = bool(params.get("filter_doublets", False))
        if filter_doublets:
            output = output.subset_obs(~predicted)

        metrics: Dict[str, Any] = {
            "n_cells_before": int(data.n_obs),
            "n_cells_after": int(output.n_obs),
            "n_doublets_marked": int(predicted.sum()),
            "n_doublets_filtered": int(predicted.sum()) if filter_doublets else 0,
            "retention_rate": float(output.n_obs / data.n_obs),
            "batch_key": batch_key,
            "batch_counts": {
                str(key): int(value)
                for key, value in (
                    data.obs[batch_key].astype(str).value_counts(dropna=False).to_dict().items()
                    if batch_key is not None
                    else {"__all__": data.n_obs}.items()
                )
            },
            "threshold_requested": threshold,
            "threshold_by_batch": threshold_by_batch,
            "threshold": float(np.mean(thresholds)) if thresholds else threshold,
            "threshold_defined": bool(thresholds) or threshold is not None,
            "n_prin_comps_requested": requested_n_pc,
            "n_prin_comps_used": effective_n_pc,
            "counts_source": counts_source,
            "counts_layer_preserved": "counts" in output.layers or counts_source == "X",
            "filter_doublets": filter_doublets,
        }

        operations = list(self.contract_operations)
        if filter_doublets:
            # Use the common contract operation name so the side-effect guard
            # can reject filtering when the caller forbids filter_cells.
            operations.append("filter_cells")

        table = pd.DataFrame(
            {
                "cell_id": data.obs["cell_id"].astype(str).to_numpy()
                if "cell_id" in data.obs
                else data.obs.index.astype(str).to_numpy(),
                "doublet_score": score,
                "predicted_doublet": predicted,
                "retained_after_filter": ~predicted if filter_doublets else np.ones(data.n_obs, dtype=bool),
            }
        )
        if batch_key is not None:
            table[batch_key] = data.obs[batch_key].astype(str).to_numpy()
        optional_uri = _optional_output_uri(contract, 1)
        optional_payload: Any = table
        if optional_uri is not None and _output_type(optional_uri) == ArtifactType.JSON:
            optional_payload = metrics
        result = _register_outputs(
            contract=contract,
            registry=registry,
            in_uri=in_uri,
            payload=output,
            capability=self,
            primary_default_path="doublet_detection/v1",
            operation="scrublet_doublet_detection",
            metrics=metrics,
            software_versions={"scanpy": _package_version("scanpy"), "scrublet_implementation": "scanpy.pp.scrublet"},
            optional_payload=optional_payload,
            optional_operation="scrublet_score_table",
        )
        # _register_outputs reads the class-level operation list.  Return the
        # actual operation list for this invocation (filtering is conditional).
        result.executed_operations = operations
        return result


# ---------------------------------------------------------------------------
# CellTypist
# ---------------------------------------------------------------------------


class CellTypistAnnotationCapability(BaseCapability):
    """Annotate cells using an explicitly supplied local CellTypist model."""

    contract_operations = [
        "validate_local_annotation_model",
        "celltypist_annotation",
        "record_annotation_conflicts",
    ]

    def __init__(self):
        super().__init__(
            capability_name="cell_annotation",
            implementation_id="celltypist_local_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA, ArtifactType.TABLE, ArtifactType.JSON],
        )
        self.legacy_aliases = {"celltypist": "celltypist_local_v1"}

    @staticmethod
    def _model_metadata(model: Any, params: Mapping[str, Any]) -> Dict[str, Any]:
        supplied = params.get("model_metadata", {})
        metadata: Dict[str, Any] = dict(supplied) if isinstance(supplied, Mapping) else {}
        model_metadata = getattr(model, "metadata", None)
        if isinstance(model_metadata, Mapping):
            for key, value in model_metadata.items():
                metadata.setdefault(str(key), _as_json_value(value))
        description = getattr(model, "description", None)
        if description is not None:
            metadata.setdefault("description", _as_json_value(description))

        aliases = {
            "species": ("model_species", "species"),
            "tissue": ("model_tissue", "tissue"),
        }
        for field, keys in aliases.items():
            if field in metadata and metadata[field] not in (None, ""):
                metadata[field] = _as_json_value(metadata[field])
                continue
            for key in keys:
                value = params.get(key)
                if value not in (None, ""):
                    metadata[field] = _as_json_value(value)
                    break
            if field not in metadata:
                for key in keys:
                    value = getattr(model, key, None)
                    if value not in (None, ""):
                        metadata[field] = _as_json_value(value)
                        break
            metadata.setdefault(field, "unknown")
        metadata["metadata_complete"] = all(
            str(metadata.get(field, "unknown")).strip().casefold() not in {"", "unknown", "none"}
            for field in ("species", "tissue")
        )
        return metadata

    @staticmethod
    def _prediction_frame(annotation: Any, n_obs: int) -> pd.DataFrame:
        labels = getattr(annotation, "predicted_labels", None)
        if labels is None and hasattr(annotation, "to_adata"):
            converted = annotation.to_adata()
            labels = getattr(converted, "obs", None)
        if labels is None:
            raise ValueError("CellTypist returned no predicted_labels table")
        if isinstance(labels, pd.Series):
            frame = labels.to_frame(name="predicted_labels")
        elif isinstance(labels, pd.DataFrame):
            frame = labels.copy(deep=True)
        else:
            array = np.asarray(labels)
            if array.ndim != 1:
                raise ValueError("CellTypist predictions are not one-dimensional")
            frame = pd.DataFrame({"predicted_labels": array})
        if len(frame) != n_obs:
            raise ValueError(
                f"CellTypist returned {len(frame)} labels for {n_obs} cells"
            )
        return frame.reset_index(drop=True)

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri, data = _load_scdata(contract, registry)
        params = contract.parameters
        model_value = params.get("model_path")
        if not model_value:
            raise ValueError(
                "celltypist_local_v1 requires an explicit local model_path; "
                "named models and automatic downloads are not allowed"
            )
        model_text = str(model_value).strip()
        if model_text.casefold().startswith(("http://", "https://")):
            raise ValueError("CellTypist model_path must be a local file, not a URL")
        model_path = Path(model_text).expanduser().resolve()
        if not model_path.is_file():
            raise FileNotFoundError(f"CellTypist local model_path does not exist: {model_path}")
        model_hash = _sha256_path(model_path)
        observed_hashes = {
            "model_path": model_hash,
            str(model_path): model_hash,
        }
        _verify_external_hashes(params, observed_hashes)

        try:
            celltypist = importlib.import_module("celltypist")
            models = importlib.import_module("celltypist.models")
        except ImportError as exc:
            raise ImportError(
                "celltypist_local_v1 requires the optional local 'celltypist' package; "
                "no model download or annotation fallback was run"
            ) from exc
        load_model = getattr(models, "Model", None)
        if load_model is None or not hasattr(load_model, "load"):
            raise ImportError("Installed CellTypist package does not expose models.Model.load")
        # CellTypist tests for '/' to distinguish a path from a downloadable
        # model name. Native Windows backslashes otherwise trigger discovery.
        model = load_model.load(model_path.as_posix())
        model_metadata = self._model_metadata(model, params)

        annotation_kwargs: Dict[str, Any] = {
            "model": model,
            "majority_voting": bool(params.get("majority_voting", False)),
            "mode": params.get("mode", "best match"),
        }
        for key in ("p_thres", "over_clustering"):
            if key in params:
                annotation_kwargs[key] = params[key]
        if "use_gpu" in params:
            annotation_kwargs["use_GPU"] = bool(params["use_gpu"])
        # Unsupported arguments must fail explicitly, never silently disappear.
        annotation = celltypist.annotate(data.to_anndata(), **annotation_kwargs)
        predictions = self._prediction_frame(annotation, data.n_obs)

        majority = bool(params.get("majority_voting", False))
        label_column = str(params.get("label_column", "majority_voting" if majority else "predicted_labels"))
        if label_column not in predictions.columns:
            if "predicted_labels" in predictions.columns:
                label_column = "predicted_labels"
            else:
                candidate_columns = [column for column in predictions.columns if predictions[column].dtype == object]
                if not candidate_columns:
                    raise ValueError("CellTypist predictions contain no usable label column")
                label_column = str(candidate_columns[0])
        candidate = predictions[label_column].astype(object).to_numpy()
        unknown = np.asarray([_label_is_unknown(value) for value in candidate], dtype=bool)

        existing_key = str(params.get("existing_cell_type_key", "cell_type"))
        annotation_key = str(params.get("annotation_key", "cell_type_celltypist"))
        if annotation_key == existing_key and not bool(params.get("use_as_cell_type", False)):
            raise ValueError(
                "annotation_key cannot equal existing_cell_type_key unless "
                "use_as_cell_type=True is explicit"
            )
        output = data.copy()
        output.obs[annotation_key] = candidate
        output.obs[f"{annotation_key}_unknown"] = unknown

        metadata_conflicts: Dict[str, Dict[str, str]] = {}
        for field in ("species", "tissue"):
            model_value_field = str(model_metadata.get(field, "unknown"))
            data_value = _single_metadata_value(data, field)
            if (
                data_value
                and model_value_field.casefold() not in {"", "unknown", "none"}
                and model_value_field.casefold() != data_value.casefold()
            ):
                metadata_conflicts[field] = {"data": data_value, "model": model_value_field}

        label_conflict = np.zeros(data.n_obs, dtype=bool)
        if existing_key in data.obs.columns:
            existing_values = data.obs[existing_key].astype(object).to_numpy()
            existing_unknown = np.asarray([_label_is_unknown(value) for value in existing_values], dtype=bool)
            label_conflict = (~existing_unknown) & (~unknown) & (
                existing_values.astype(str) != candidate.astype(str)
            )
            output.obs["cell_type_existing"] = existing_values
            output.obs[f"{annotation_key}_conflict"] = label_conflict
        else:
            existing_values = np.full(data.n_obs, None, dtype=object)
            output.obs[f"{annotation_key}_conflict"] = False

        use_as_cell_type = bool(params.get("use_as_cell_type", False))
        if use_as_cell_type:
            if existing_key in data.obs.columns:
                # Unknown model labels, metadata conflicts, and row-level label
                # conflicts retain the pre-existing annotation.
                final_values = existing_values.copy()
                eligible = (~unknown) & (~label_conflict) & (not bool(metadata_conflicts))
                final_values[eligible] = candidate[eligible]
                output.obs[existing_key] = final_values
            elif not metadata_conflicts:
                output.obs[existing_key] = candidate

        output.uns["celltypist"] = {
            "implementation_id": self.implementation_id,
            "model_path": str(model_path),
            "model_sha256": model_hash,
            "model_metadata": model_metadata,
            "annotation_key": annotation_key,
            "existing_cell_type_key": existing_key if existing_key in data.obs.columns else None,
            "use_as_cell_type": use_as_cell_type,
            "metadata_conflicts": metadata_conflicts,
            "unknown_label_count": int(unknown.sum()),
            "label_conflict_count": int(label_conflict.sum()),
        }

        metrics: Dict[str, Any] = {
            "n_cells_before": int(data.n_obs),
            "n_cells_after": int(output.n_obs),
            "n_annotations": int((~unknown).sum()),
            "n_unknown_annotations": int(unknown.sum()),
            "n_label_conflicts": int(label_conflict.sum()),
            "metadata_conflicts": metadata_conflicts,
            "model_path": str(model_path),
            "model_sha256": model_hash,
            "model_metadata": model_metadata,
            "annotation_key": annotation_key,
            "use_as_cell_type": use_as_cell_type,
            "existing_cell_type_preserved": existing_key in data.obs.columns,
        }
        table = pd.DataFrame({
            "cell_id": data.obs["cell_id"].astype(str).to_numpy()
            if "cell_id" in data.obs
            else data.obs.index.astype(str).to_numpy(),
            "celltypist_label": candidate,
            "celltypist_unknown": unknown,
            "celltypist_conflict": label_conflict,
        })
        for column in predictions.columns:
            if column == label_column:
                continue
            values = predictions[column].to_numpy()
            if np.asarray(values).ndim == 1 and len(values) == data.n_obs:
                table[f"celltypist_{column}"] = values
        optional_uri = _optional_output_uri(contract, 1)
        optional_payload = table
        if optional_uri is not None and _output_type(optional_uri) == ArtifactType.JSON:
            optional_payload = metrics
        return _register_outputs(
            contract=contract,
            registry=registry,
            in_uri=in_uri,
            payload=output,
            capability=self,
            primary_default_path="cell_annotation/v1",
            operation="celltypist_annotation",
            metrics=metrics,
            software_versions={"celltypist": _package_version("celltypist")},
            optional_payload=optional_payload,
            optional_operation="celltypist_annotation_table",
        )


# ---------------------------------------------------------------------------
# CellBender CLI adapter
# ---------------------------------------------------------------------------


class CellBenderBackgroundRemovalCapability(BaseCapability):
    """Run an explicit CellBender executable and retain corrected counts."""

    contract_operations = [
        "validate_unfiltered_droplet_input",
        "run_cellbender_cli",
        "validate_cell_gene_alignment",
        "store_corrected_counts_layer",
    ]

    def __init__(self):
        super().__init__(
            capability_name="background_removal",
            implementation_id="cellbender_cli_v1",
            implementation_type=ImplementationType.CONTAINER,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA, ArtifactType.JSON],
        )
        self.legacy_aliases = {"cellbender": "cellbender_cli_v1"}

    @staticmethod
    def _resolve_executable(value: Any) -> Tuple[str, str]:
        if not value:
            raise ValueError(
                "cellbender_cli_v1 requires an explicit executable parameter; "
                "a missing executable cannot be replaced with a synthetic correction"
            )
        raw = str(value)
        resolved = Path(raw).expanduser().resolve()
        if resolved.is_file():
            return str(resolved), _sha256_path(resolved)
        located = shutil.which(raw)
        if located is None:
            raise FileNotFoundError(f"CellBender executable is not available: {raw}")
        located_path = Path(located).resolve()
        return str(located_path), _sha256_path(located_path) if located_path.is_file() else "unavailable"

    @staticmethod
    def _read_output(path: Path) -> Any:
        try:
            import anndata as ad
        except ImportError as exc:
            raise ImportError(
                "Reading CellBender output requires the optional 'anndata' package"
            ) from exc
        suffix = path.suffix.casefold()
        suffixes = "".join(path.suffixes).casefold()
        if suffix == ".h5ad" or suffixes.endswith(".h5ad.gz"):
            return ad.read_h5ad(path)
        try:
            scanpy = importlib.import_module("scanpy")
        except ImportError as exc:
            raise ImportError(
                "CellBender .h5 or 10x directory output requires Scanpy's reader"
            ) from exc
        if path.is_dir():
            return scanpy.read_10x_mtx(path, var_names="gene_symbols", make_unique=False)
        if suffix == ".h5":
            return scanpy.read_10x_h5(path, genome=None)
        raise ValueError(
            "Unsupported CellBender output format; use .h5ad, .h5, or a 10x matrix directory"
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri, data = _load_scdata(contract, registry)
        params = contract.parameters
        raw_value = params.get("unfiltered_input_path")
        if not raw_value:
            raise ValueError(
                "cellbender_cli_v1 requires an explicit unfiltered_input_path; "
                "filtered artifact X cannot stand in for CellBender's droplet input"
            )
        raw_path = Path(str(raw_value)).expanduser().resolve()
        if not raw_path.exists():
            raise FileNotFoundError(f"CellBender unfiltered_input_path does not exist: {raw_path}")
        raw_hash = _sha256_path(raw_path)
        executable, executable_hash = self._resolve_executable(params.get("executable"))
        observed_hashes = {
            "unfiltered_input_path": raw_hash,
            str(raw_path): raw_hash,
            "executable": executable_hash,
            str(Path(executable).resolve()): executable_hash,
        }

        output_value = params.get("output_path")
        if not output_value:
            raise ValueError(
                "cellbender_cli_v1 requires an explicit output_path so the external result "
                "and its run path can be audited"
            )
        output_path = Path(str(output_value)).expanduser().resolve()
        if output_path == raw_path:
            raise ValueError("CellBender output_path must differ from unfiltered_input_path")
        if output_path.suffix.casefold() != ".h5":
            raise ValueError(
                "CellBender remove-background output_path must end with .h5 "
                "(the CellRanger v3 output format)"
            )
        if output_path.exists() and not bool(params.get("allow_overwrite", False)):
            raise FileExistsError(
                f"CellBender output_path already exists: {output_path}; "
                "set allow_overwrite=True explicitly to reuse it"
            )
        if not output_path.parent.exists():
            raise FileNotFoundError(f"CellBender output parent directory does not exist: {output_path.parent}")

        run_cwd = Path(str(params.get("run_cwd", Path.cwd()))).expanduser().resolve()
        if not run_cwd.is_dir():
            raise FileNotFoundError(f"CellBender run_cwd does not exist: {run_cwd}")
        subcommand = str(params.get("subcommand", "remove-background"))
        if subcommand != "remove-background":
            raise ValueError("cellbender_cli_v1 only supports the remove-background subcommand")
        command: List[str] = [
            executable,
            subcommand,
            "--input",
            str(raw_path),
            "--output",
            str(output_path),
        ]
        extra_args = params.get("extra_args", params.get("cli_args", []))
        if extra_args is None:
            extra_args = []
        if not isinstance(extra_args, (list, tuple)) or any(not isinstance(item, str) for item in extra_args):
            raise ValueError("CellBender extra_args/cli_args must be a list of strings")
        if any(
            item in {"--input", "--output"}
            or item.startswith("--input=")
            or item.startswith("--output=")
            for item in extra_args
        ):
            raise ValueError(
                "extra_args/cli_args cannot override the adapter's explicit --input or --output"
            )
        checkpoint_resources = _cellbender_checkpoint_resources(extra_args, run_cwd)
        for key, checkpoint in checkpoint_resources:
            checkpoint_hash = _sha256_path(checkpoint)
            observed_hashes[key] = checkpoint_hash
            observed_hashes[str(checkpoint.resolve())] = checkpoint_hash
        _verify_external_hashes(params, observed_hashes)
        command.extend(str(item) for item in extra_args)
        timeout = params.get("timeout_sec", 3600)
        try:
            completed = subprocess.run(
                command,
                cwd=str(run_cwd),
                capture_output=True,
                text=True,
                check=False,
                timeout=float(timeout),
            )
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"CellBender executable failed to launch: {executable}") from exc
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"CellBender CLI timed out after {timeout} seconds in {run_cwd}"
            ) from exc
        stdout = _clean_subprocess_text(getattr(completed, "stdout", ""))
        stderr = _clean_subprocess_text(getattr(completed, "stderr", ""))
        if int(completed.returncode) != 0:
            raise RuntimeError(
                f"CellBender CLI failed with return code {completed.returncode}; "
                f"cwd={run_cwd}; stderr={stderr[-2000:]}"
            )
        if not output_path.exists():
            raise FileNotFoundError(
                f"CellBender exited successfully but produced no output at {output_path}"
            )
        output_hash = _sha256_path(output_path)
        # CellBender may exit zero while report generation is unavailable or
        # while its report flags convergence concerns.  Capture those
        # sidecars before publishing the AnnData artifact so the caveat is
        # carried by both the artifact and the optional run report.
        sidecars = _collect_cellbender_sidecars(
            output_path,
            params,
            stdout=stdout,
            stderr=stderr,
        )

        external = self._read_output(output_path)
        if not hasattr(external, "X") or not hasattr(external, "obs") or not hasattr(external, "var"):
            raise ValueError("CellBender output is not an AnnData-like matrix")
        input_cells, input_genes = _cell_and_gene_names(data, params)
        output_cells, output_genes = _external_anndata_names(external, params)
        corrected, mapping = _map_external_matrix_to_target(
            external.X,
            output_cells,
            output_genes,
            input_cells,
            input_genes,
        )
        values = _matrix_values(corrected)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError("CellBender corrected output contains nonfinite or negative counts")
        if corrected.shape != data.shape:
            raise ValueError("CellBender corrected output shape changed during alignment")

        output = data.copy()
        if "counts" not in output.layers:
            output.layers["counts"] = _copy_matrix(output.X)
        original_counts = _copy_matrix(output.layers["counts"])
        output.layers["corrected_counts"] = corrected
        if original_counts.shape != output.layers["counts"].shape:
            raise ValueError("Original counts layer changed shape while storing CellBender output")

        report: Dict[str, Any] = {
            "implementation_id": self.implementation_id,
            "executable": executable,
            "executable_sha256": executable_hash,
            "unfiltered_input_path": str(raw_path),
            "unfiltered_input_sha256": raw_hash,
            "output_path": str(output_path),
            "output_sha256": output_hash,
            "matrix_source": "cellbender_cli_output",
            "run_cwd": str(run_cwd),
            "command": command,
            "returncode": int(completed.returncode),
            "stdout": stdout[-4000:],
            "stderr": stderr[-4000:],
            "cell_gene_alignment": mapping,
            "raw_counts_layer": "counts",
            "corrected_counts_layer": "corrected_counts",
            "original_counts_preserved": True,
            "checkpoint_paths": [str(path) for _, path in checkpoint_resources],
            "checkpoint_sha256": {
                key: observed_hashes[key] for key, _ in checkpoint_resources
            },
        }
        report.update(sidecars)
        output.uns["cellbender"] = report
        metrics: Dict[str, Any] = {
            "n_cells_before": int(data.n_obs),
            "n_cells_after": int(output.n_obs),
            "n_genes": int(output.n_vars),
            "unfiltered_input_path": str(raw_path),
            "unfiltered_input_sha256": raw_hash,
            "output_path": str(output_path),
            "output_sha256": output_hash,
            "matrix_source": "cellbender_cli_output",
            "executable": executable,
            "executable_sha256": executable_hash,
            "run_cwd": str(run_cwd),
            "returncode": int(completed.returncode),
            "cells_aligned": True,
            "genes_aligned": True,
            "cell_gene_mapping": mapping,
            "original_counts_preserved": True,
            "corrected_counts_layer": "corrected_counts",
            "checkpoint_paths": [str(path) for _, path in checkpoint_resources],
            "checkpoint_sha256": {
                key: observed_hashes[key] for key, _ in checkpoint_resources
            },
        }
        metrics.update(sidecars)
        optional_uri = _optional_output_uri(contract, 1)
        optional_payload = metrics
        return _register_outputs(
            contract=contract,
            registry=registry,
            in_uri=in_uri,
            payload=output,
            capability=self,
            primary_default_path="background_removal/v1",
            operation="cellbender_background_removal",
            metrics=metrics,
            software_versions={"cellbender_executable": executable},
            optional_payload=optional_payload,
            optional_operation="cellbender_run_report",
        )


# Friendly aliases make the implementation discoverable for callers that use
# the method name instead of the longer capability class name.
ScanpyScrubletCapability = ScrubletDoubletCapability
DoubletDetectionCapability = ScrubletDoubletCapability
CellTypistCapability = CellTypistAnnotationCapability
CellTypeAnnotationCapability = CellTypistAnnotationCapability
CellBenderCapability = CellBenderBackgroundRemovalCapability
BackgroundRemovalCapability = CellBenderBackgroundRemovalCapability


__all__ = [
    "ScrubletDoubletCapability",
    "ScanpyScrubletCapability",
    "DoubletDetectionCapability",
    "CellTypistAnnotationCapability",
    "CellTypistCapability",
    "CellTypeAnnotationCapability",
    "CellBenderBackgroundRemovalCapability",
    "CellBenderCapability",
    "BackgroundRemovalCapability",
]
