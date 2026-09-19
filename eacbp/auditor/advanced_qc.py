"""Validate cell identity, count preservation and explicit QC decisions."""
import hashlib
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity
from eacbp.capabilities.sc_data import SCData


def same_matrix(a, b):
    if a.shape != b.shape:
        return False
    if sparse.issparse(a) or sparse.issparse(b):
        difference = sparse.csr_matrix(a) - sparse.csr_matrix(b)
        return difference.nnz == 0
    return np.array_equal(a, b)


def _sha256_path(path):
    """Hash the exact external file or deterministic directory tree."""

    path = Path(path).expanduser().resolve()
    digest = hashlib.sha256()
    if path.is_file():
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    if path.is_dir():
        files = sorted(item for item in path.rglob("*") if item.is_file())
        if not files:
            raise ValueError(f"External directory is empty: {path}")
        for item in files:
            relative = item.relative_to(path).as_posix().encode("utf-8")
            digest.update(len(relative).to_bytes(8, "big"))
            digest.update(relative)
            with item.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
        return digest.hexdigest()
    raise FileNotFoundError(f"External resource does not exist: {path}")


def _hash_value(value):
    return str(value).strip().lower().removeprefix("sha256:")


def _cellbender_checkpoint_resources(parameters):
    """Resolve the adapter's local ``--checkpoint`` arguments for auditing."""

    if "extra_args" in parameters:
        raw_args = parameters.get("extra_args")
    elif "cli_args" in parameters:
        raw_args = parameters.get("cli_args")
    else:
        return []
    if not isinstance(raw_args, (list, tuple)) or any(not isinstance(item, str) for item in raw_args):
        raise ValueError("CellBender extra_args/cli_args must be a list of strings")
    run_cwd = Path(str(parameters.get("run_cwd", Path.cwd()))).expanduser().resolve()
    if not run_cwd.is_dir():
        raise FileNotFoundError(f"CellBender run_cwd does not exist: {run_cwd}")

    resources = []
    index = 0
    while index < len(raw_args):
        argument = raw_args[index]
        checkpoint_value = None
        if argument == "--checkpoint":
            if index + 1 >= len(raw_args):
                raise ValueError("CellBender extra_args/cli_args has --checkpoint without a path")
            checkpoint_value = raw_args[index + 1]
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


def _warning_strings(value):
    """Flatten quality warnings that may round-trip through AnnData arrays."""

    if value is None:
        return []
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return _warning_strings(value.item())
        flattened = []
        for item in value:
            flattened.extend(_warning_strings(item))
        return flattened
    if isinstance(value, np.generic):
        return _warning_strings(value.item())
    if isinstance(value, (list, tuple, set)):
        flattened = []
        for item in value:
            flattened.extend(_warning_strings(item))
        return flattened
    text = str(value).strip()
    return [text] if text else []


def _matrix_ids(adata, params):
    cell_key = str(params.get("cell_id_key", "cell_id"))
    gene_key = str(params.get("gene_id_key", "gene_name"))
    cells = (
        adata.obs[cell_key].astype(str).to_numpy()
        if cell_key in adata.obs.columns
        else np.asarray(getattr(adata, "obs_names", adata.obs.index).astype(str))
    )
    genes = (
        adata.var[gene_key].astype(str).to_numpy()
        if gene_key in adata.var.columns
        else np.asarray(getattr(adata, "var_names", adata.var.index).astype(str))
    )
    if len(cells) != len(set(cells)):
        raise ValueError("External matrix has duplicate cell identifiers")
    if len(genes) != len(set(genes)):
        raise ValueError("External matrix has duplicate gene identifiers")
    return cells, genes


def _map_external_matrix(matrix, external_cells, external_genes, target_cells, target_genes):
    """Map a full CellBender output onto target IDs, rejecting omissions."""

    external_cells = np.asarray([str(value) for value in external_cells], dtype=object)
    external_genes = np.asarray([str(value) for value in external_genes], dtype=object)
    target_cells = np.asarray([str(value) for value in target_cells], dtype=object)
    target_genes = np.asarray([str(value) for value in target_genes], dtype=object)
    if len(external_cells) != len(set(external_cells)):
        raise ValueError("CellBender output has duplicate cell identifiers")
    if len(external_genes) != len(set(external_genes)):
        raise ValueError("CellBender output has duplicate gene identifiers")
    cell_index = {value: index for index, value in enumerate(external_cells)}
    gene_index = {value: index for index, value in enumerate(external_genes)}
    missing_cells = [value for value in target_cells if value not in cell_index]
    missing_genes = [value for value in target_genes if value not in gene_index]
    if missing_cells:
        raise ValueError(f"CellBender output is missing target cell IDs: {missing_cells[:10]}")
    if missing_genes:
        raise ValueError(f"CellBender output is missing target gene IDs: {missing_genes[:10]}")
    rows = np.asarray([cell_index[value] for value in target_cells], dtype=int)
    columns = np.asarray([gene_index[value] for value in target_genes], dtype=int)
    if sparse.issparse(matrix):
        mapped = matrix[rows][:, columns].copy()
    else:
        mapped = np.asarray(matrix)[np.ix_(rows, columns)].copy()
    return mapped, {
        "target_n_cells": int(len(target_cells)),
        "target_n_genes": int(len(target_genes)),
        "external_n_cells": int(len(external_cells)),
        "external_n_genes": int(len(external_genes)),
        "mapped_n_cells": int(len(rows)),
        "mapped_n_genes": int(len(columns)),
        "extra_external_cells": int(len(external_cells) - len(target_cells)),
        "extra_external_genes": int(len(external_genes) - len(target_genes)),
        "cell_order_exact": bool(np.array_equal(external_cells, target_cells)),
        "gene_order_exact": bool(np.array_equal(external_genes, target_genes)),
        "cell_ids_unique": True,
        "gene_ids_unique": True,
    }


def _read_cellbender_output(path):
    """Read the external CellBender matrix without trusting capability state."""

    path = Path(path)
    if path.suffix.casefold() == ".h5":
        import scanpy as sc

        return sc.read_10x_h5(path, genome=None)
    if path.suffix.casefold() == ".h5ad":
        import anndata as ad

        return ad.read_h5ad(path)
    if path.is_dir():
        import scanpy as sc

        return sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=False)
    raise ValueError("Unsupported CellBender output format for independent audit")


class AdvancedQCValidator(BaseAuditor):
    def __init__(self):
        super().__init__("advanced_qc_validator")

    def audit(self, contract, result, registry):
        report = ValidationReport(auditor_name=self.auditor_name, target_task_id=contract.task_id)
        if contract.capability not in {"doublet_detection", "cell_annotation", "background_removal"}:
            return report
        errors = []
        cellbender_quality_info = {}
        try:
            _, source = registry.get(contract.input_artifacts[0])
            _, output = registry.get(result.output_artifacts[0])
            source = source if isinstance(source, SCData) else SCData.from_dict(source)
            output = output if isinstance(output, SCData) else SCData.from_dict(output)
            params = contract.parameters
            if not source.obs.index.is_unique or not output.obs.index.is_unique:
                raise ValueError("Unique observation IDs are required to audit cell correspondence")
            positions = source.obs.index.get_indexer(output.obs.index)
            if (positions < 0).any() or not source.var.equals(output.var):
                raise ValueError("QC output has unknown cells or altered gene identity")
            filtering = contract.capability == "doublet_detection" and bool(params.get("filter_doublets", False))
            if not filtering and not source.obs.index.equals(output.obs.index):
                errors.append("Cell identities/order changed without explicit filtering")
            if not same_matrix(source.X[positions], output.X):
                errors.append("QC/annotation unexpectedly changed expression values")
            for layer, matrix in source.layers.items():
                if layer not in output.layers or not same_matrix(matrix[positions], output.layers[layer]):
                    errors.append(f"Original layer changed: {layer}")
            if "counts" in output.layers:
                raw = source.layers.get("counts", source.X)
                if not same_matrix(raw[positions], output.layers["counts"]):
                    errors.append("Raw counts changed")
            if contract.capability == "doublet_detection":
                if not {"doublet_score", "predicted_doublet"}.issubset(output.obs):
                    raise ValueError("Doublet outputs are absent")
                scores = pd.to_numeric(output.obs.doublet_score, errors="coerce")
                if not scores.between(0, 1).all() or not output.obs.predicted_doublet.isin([True, False]).all():
                    errors.append("Invalid doublet scores or labels")
                info = output.uns.get("eacbp_scrublet", {})
                thresholds = info.get("threshold_by_batch", {})
                batch = info.get("batch_key")
                labels = output.obs[batch].astype(str) if batch else pd.Series("__all__", index=output.obs.index)
                for label in labels.unique():
                    threshold = thresholds.get(str(label), info.get("threshold_requested"))
                    if threshold is None or not np.isfinite(float(threshold)):
                        errors.append(f"Unresolved doublet threshold for batch {label}")
                        continue
                    rows = labels == label
                    if not np.array_equal(scores.loc[rows] > float(threshold), output.obs.loc[rows, "predicted_doublet"].astype(bool)):
                        errors.append("Doublet labels do not match recorded threshold")
                if filtering and output.obs.predicted_doublet.any():
                    errors.append("Filtered output retains predicted doublets")
                if filtering:
                    tables = [registry.get(uri)[1] for uri in result.output_artifacts[1:]
                              if uri.startswith("table://")]
                    if not tables:
                        errors.append("Filtering requires a complete pre-filter score table")
                    else:
                        table = tables[0]
                        ids = (source.obs["cell_id"].astype(str) if "cell_id" in source.obs
                               else pd.Series(source.obs.index.astype(str), index=source.obs.index))
                        if not ids.is_unique or table.cell_id.astype(str).tolist() != ids.tolist():
                            raise ValueError("Pre-filter scores must cover every original cell in order")
                        all_scores = pd.to_numeric(table.doublet_score, errors="coerce")
                        if not all_scores.between(0, 1).all() or not table.predicted_doublet.isin([True, False]).all():
                            raise ValueError("Invalid pre-filter doublet scores or labels")
                        keep = ~table.predicted_doublet.to_numpy(dtype=bool)
                        if not np.array_equal(table.retained_after_filter.to_numpy(), keep):
                            errors.append("Recorded retention differs from doublet labels")
                        if not source.obs.index[keep].equals(output.obs.index):
                            errors.append("Filtered cells differ from the complete score table")
                        all_labels = source.obs[batch].astype(str).to_numpy() if batch else np.repeat("__all__", source.n_obs)
                        for label in np.unique(all_labels):
                            threshold = thresholds.get(str(label), info.get("threshold_requested"))
                            rows = all_labels == label
                            if threshold is None or not np.array_equal(
                                all_scores.to_numpy()[rows] > float(threshold),
                                table.predicted_doublet.to_numpy(dtype=bool)[rows],
                            ):
                                errors.append("Pre-filter labels do not match recorded thresholds")
            elif contract.capability == "cell_annotation":
                info = output.uns.get("celltypist", {})
                model_path = Path(params["model_path"])
                actual = hashlib.sha256(model_path.read_bytes()).hexdigest()
                if info.get("model_sha256") != actual:
                    errors.append("Model hash differs from recorded annotation resource")
                key = params.get("annotation_key", "cell_type_celltypist")
                if key not in output.obs or f"{key}_unknown" not in output.obs:
                    errors.append("Annotation or uncertainty labels are absent")
                existing = params.get("existing_cell_type_key", "cell_type")
                if not params.get("use_as_cell_type", False) and existing in source.obs and not source.obs[existing].equals(output.obs[existing]):
                    errors.append("Existing annotation changed without explicit request")
            else:
                corrected = output.layers.get("corrected_counts")
                if corrected is None:
                    raise ValueError("Corrected counts layer is absent")
                values = corrected.data if sparse.issparse(corrected) else np.asarray(corrected)
                if not np.isfinite(values).all() or (values < 0).any():
                    errors.append("Corrected counts are invalid")
                info = output.uns.get("cellbender", {})
                if not isinstance(info, dict):
                    errors.append("CellBender provenance record is absent or malformed")
                else:
                    cellbender_quality_info = info
                    required = {
                        "unfiltered_input_path",
                        "unfiltered_input_sha256",
                        "executable",
                        "executable_sha256",
                        "output_path",
                        "output_sha256",
                        "run_cwd",
                        "matrix_source",
                        "cell_gene_alignment",
                    }
                    missing = sorted(required - set(info))
                    if missing:
                        errors.append(f"CellBender provenance fields are missing: {missing}")
                    if info.get("matrix_source") != "cellbender_cli_output":
                        errors.append("Corrected counts are not identified as CellBender CLI output")
                    if info.get("corrected_counts_layer") != "corrected_counts":
                        errors.append("CellBender corrected layer provenance is invalid")

                    path_values = {
                        "unfiltered_input_path": info.get("unfiltered_input_path"),
                        "executable": info.get("executable"),
                        "output_path": info.get("output_path"),
                        "run_cwd": info.get("run_cwd"),
                    }
                    expected_paths = {
                        key: params.get(key)
                        for key in ("unfiltered_input_path", "executable", "output_path", "run_cwd")
                        if params.get(key) is not None
                    }
                    for key, expected in expected_paths.items():
                        if key == "executable" and not Path(str(expected)).exists():
                            # The adapter resolves PATH commands to an absolute
                            # executable path; resolve the expected command only
                            # for comparison when it is still available.
                            import shutil

                            resolved = shutil.which(str(expected))
                            expected = resolved or expected
                        try:
                            expected_resolved = str(Path(str(expected)).expanduser().resolve())
                            actual_resolved = str(Path(str(path_values.get(key))).expanduser().resolve())
                            if expected_resolved != actual_resolved:
                                errors.append(f"CellBender {key} differs from the task contract")
                        except (TypeError, ValueError):
                            errors.append(f"CellBender {key} provenance is not a valid path")

                    observed_hashes = {}
                    for key, value in path_values.items():
                        if value is None or key == "run_cwd":
                            continue
                        try:
                            observed_hashes[key] = _sha256_path(value)
                        except (OSError, ValueError) as exc:
                            errors.append(f"CellBender {key} cannot be hashed: {exc}")
                    try:
                        checkpoint_resources = _cellbender_checkpoint_resources(params)
                    except (OSError, TypeError, ValueError) as exc:
                        checkpoint_resources = []
                        errors.append(f"CellBender checkpoint cannot be hashed: {exc}")
                    for key, checkpoint in checkpoint_resources:
                        try:
                            checkpoint_hash = _sha256_path(checkpoint)
                        except (OSError, ValueError) as exc:
                            errors.append(f"CellBender {key} cannot be hashed: {exc}")
                        else:
                            observed_hashes[key] = checkpoint_hash
                            observed_hashes[str(checkpoint.resolve())] = checkpoint_hash
                    recorded_hashes = {
                        "unfiltered_input_path": info.get("unfiltered_input_sha256"),
                        "executable": info.get("executable_sha256"),
                        "output_path": info.get("output_sha256"),
                    }
                    for key, recorded in recorded_hashes.items():
                        actual = observed_hashes.get(key)
                        if actual is None:
                            continue
                        if recorded is None or _hash_value(recorded) != _hash_value(actual):
                            errors.append(f"CellBender {key} hash differs from provenance")

                    # Sidecar files are independent evidence emitted by the
                    # external process.  A missing report is a quality caveat,
                    # while a changed sidecar with a recorded hash is an
                    # integrity failure just like a changed output matrix.
                    for kind in ("report", "metrics", "log"):
                        sidecar_path = info.get(f"{kind}_path")
                        if not sidecar_path:
                            continue
                        sidecar_file = Path(str(sidecar_path)).expanduser().resolve()
                        sidecar_hash = info.get(f"{kind}_sha256")
                        sidecar_available = bool(info.get(f"{kind}_available"))
                        if sidecar_file.is_file():
                            try:
                                actual_sidecar_hash = _sha256_path(sidecar_file)
                            except (OSError, ValueError) as exc:
                                errors.append(f"CellBender {kind} sidecar cannot be hashed: {exc}")
                            else:
                                if sidecar_hash and _hash_value(sidecar_hash) != _hash_value(actual_sidecar_hash):
                                    errors.append(f"CellBender {kind} sidecar hash differs from provenance")
                                if not sidecar_available:
                                    errors.append(f"CellBender {kind} sidecar availability is inconsistent")
                        elif sidecar_available:
                            errors.append(f"CellBender {kind} sidecar is recorded as available but is missing")

                    expected_hashes = params.get("external_resource_sha256", {})
                    if expected_hashes is not None:
                        if not isinstance(expected_hashes, dict):
                            errors.append("external_resource_sha256 is not a mapping")
                        else:
                            for key, expected in expected_hashes.items():
                                actual = observed_hashes.get(str(key))
                                if actual is None:
                                    path_key = str(Path(str(key)).expanduser().resolve())
                                    for candidate_key, candidate_value in path_values.items():
                                        if candidate_value is not None and str(Path(str(candidate_value)).expanduser().resolve()) == path_key:
                                            actual = observed_hashes.get(candidate_key)
                                            break
                                if actual is None or _hash_value(actual) != _hash_value(expected):
                                    errors.append(f"CellBender external hash mismatch for {key}")

                    alignment = info.get("cell_gene_alignment", {})
                    if not isinstance(alignment, dict):
                        errors.append("CellBender cell/gene alignment record is malformed")
                    else:
                        if alignment.get("mapped_n_cells") != output.n_obs or alignment.get("mapped_n_genes") != output.n_vars:
                            errors.append("CellBender alignment record does not match target artifact dimensions")

                    # Hash and path checks above are independent of the
                    # capability's in-memory object.  If the external output
                    # is still available, read it again and verify that the
                    # corrected layer is exactly the ID-mapped external matrix.
                    output_path = info.get("output_path")
                    if output_path and Path(str(output_path)).exists():
                        try:
                            external = _read_cellbender_output(output_path)
                            source_cells, source_genes = _matrix_ids(source, params)
                            external_cells, external_genes = _matrix_ids(external, params)
                            mapped, mapping = _map_external_matrix(
                                external.X,
                                external_cells,
                                external_genes,
                                source_cells,
                                source_genes,
                            )
                            if not same_matrix(mapped, corrected):
                                errors.append("corrected_counts does not match the ID-mapped CellBender output")
                            recorded_mapping = alignment if isinstance(alignment, dict) else {}
                            for key, value in mapping.items():
                                if recorded_mapping.get(key) != value:
                                    errors.append(f"CellBender alignment provenance differs for {key}")
                        except Exception as exc:
                            errors.append(f"CellBender external matrix audit failed: {type(exc).__name__}: {exc}")
                    else:
                        errors.append("CellBender output path is unavailable for independent matrix audit")
        except Exception as exc:
            errors.append(f"{type(exc).__name__}: {exc}")
        report.add_check("advanced_qc_integrity", not errors, ValidationSeverity.ERROR,
                         "; ".join(errors) if errors else "Checked cell/gene correspondence, preserved matrices and explicit QC/annotation outputs.")
        if contract.capability == "background_removal":
            quality_warnings = _warning_strings(cellbender_quality_info.get("quality_warnings", []))
            report_available = bool(cellbender_quality_info.get("report_available"))
            evidence_recorded = any(
                bool(cellbender_quality_info.get(f"{kind}_available"))
                for kind in ("report", "metrics", "log")
            )
            quality_evidence_ok = bool(report_available and not quality_warnings)
            if quality_evidence_ok:
                quality_message = (
                    "CellBender report was retained without recorded warnings; "
                    "this is execution evidence only and does not certify scientific quality."
                )
            elif quality_warnings:
                quality_message = (
                    "CellBender execution completed, but quality evidence contains warnings; "
                    "finite corrected counts and exit code 0 are not a scientific quality certification. "
                    + "; ".join(str(item) for item in quality_warnings[:8])
                )
            elif not evidence_recorded:
                quality_message = (
                    "CellBender report/metrics/log evidence was not retained; "
                    "finite corrected counts and exit code 0 are not sufficient for scientific quality."
                )
            else:
                quality_message = (
                    "CellBender sidecar evidence is incomplete; inspect the external report before use."
                )
            report.add_check(
                "cellbender_quality_evidence",
                quality_evidence_ok,
                ValidationSeverity.WARNING,
                quality_message,
                metrics={
                    "report_available": report_available,
                    "evidence_recorded": evidence_recorded,
                    "quality_status": cellbender_quality_info.get("quality_status"),
                    "scientific_quality_certified": False,
                    "warning_count": len(quality_warnings),
                },
            )
        return report
