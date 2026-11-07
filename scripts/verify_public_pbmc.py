"""Run an auditable PBMC3k clustering check against a processed reference.

The raw H5AD is never modified. A lossless working copy receives an explicit
``layers['counts']`` layer and a gene-symbol column for the standard pipeline.
The processed reference labels are read only after EACBP has completed, and
are used only to compute barcode-aligned clustering metrics. Both public
input files must match the checksum-pinned download catalog.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePath, PureWindowsPath
import re
import sys
from typing import Any, Mapping
from uuid import uuid4

import numpy as np
import pandas as pd
from scipy import sparse

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

if __package__:
    from .fetch_public_validation_data import DATASETS
else:
    from fetch_public_validation_data import DATASETS

PUBLIC_INPUT_SHA256 = {spec.key: spec.sha256 for spec in DATASETS}


# Fixed before the run; these markers do not use reference labels or cluster
# assignments. They provide coarse PBMC labels for the standard marker scorer.
PBMC_MARKERS = {
    "B_cells": ["MS4A1", "CD79A", "CD79B", "CD74"],
    "CD4_T_cells": ["IL7R", "LTB", "CCR7", "CD3D", "CD3E"],
    "CD8_T_cells": ["CD8A", "CD8B", "CCL5", "GZMK", "CD3D"],
    "NK_cells": ["NKG7", "GNLY", "PRF1", "KLRD1", "GZMB"],
    "CD14_monocytes": ["LYZ", "S100A8", "S100A9", "FCN1", "CTSS"],
    "FCGR3A_monocytes": ["FCGR3A", "MS4A7", "LST1", "IFITM3"],
    "Dendritic_cells": ["FCER1A", "CD1C", "CLEC10A", "CST3"],
    "Platelets": ["PPBP", "PF4", "SDPR"],
}
REFERENCE_LABEL_KEY = "louvain"

# Explicit human PBMC aliases used only after clustering for label evaluation.
# Labels outside this map stay unmapped; no nearest-class or expression-based
# relabeling is performed from the reference.
PBMC_CELL_TYPE_LABEL_MAP = {
    "CD4 T cells": ("CD4 T cells", "CD4+ T cells", "CD4 T lymphocytes", "CD4_T_cells"),
    "CD8 T cells": ("CD8 T cells", "CD8+ T cells", "CD8 T lymphocytes", "CD8_T_cells"),
    "B cells": ("B cells", "B cell", "B lymphocytes", "B_cells"),
    "NK cells": ("NK cells", "natural killer cells", "natural killer lymphocytes", "NK_cells"),
    "CD14+ Monocytes": ("CD14+ Monocytes", "CD14 Monocytes", "CD14+ Mono", "CD14_monocytes"),
    "FCGR3A+ Monocytes": (
        "FCGR3A+ Monocytes", "FCGR3A Monocytes", "CD16+ Monocytes", "non-classical Monocytes",
        "FCGR3A_monocytes",
    ),
    "Dendritic cells": ("Dendritic cells", "Dendritic cell", "DCs", "Dendritic_cells"),
    "Megakaryocytes": ("Megakaryocytes", "Megakaryocyte", "Platelets", "Platelet cells"),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix_values(matrix: Any) -> np.ndarray:
    return matrix.data if sparse.issparse(matrix) else np.asarray(matrix)


def _validate_integer_counts(matrix: Any, label: str) -> None:
    values = _matrix_values(matrix)
    if values.ndim != 1 and values.ndim != 2:
        raise ValueError(f"{label} must be a two-dimensional count matrix")
    if values.size == 0:
        raise ValueError(f"{label} is empty")
    if values.dtype.kind not in "biuf" or not np.isfinite(values).all():
        raise ValueError(f"{label} must contain finite numeric counts")
    if np.any(values < 0) or (values.dtype.kind == "f" and not np.equal(values, np.floor(values)).all()):
        raise ValueError(f"{label} must contain nonnegative integer-valued raw counts; scaled/log data are not accepted")


def _gene_name_column(adata: Any) -> tuple[str, str]:
    """Ensure the exact column consumed by EACBP's marker scorer is present."""
    columns = list(adata.var.columns)
    candidates = ("gene_name", "gene_symbols", "gene_symbol", "symbol", "feature_name", "gene_short_name")
    for column in candidates:
        if column not in columns:
            continue
        values = adata.var[column]
        if values.isna().any():
            continue
        labels = values.astype(str).str.strip()
        if labels.empty or labels.eq("").any() or labels.str.casefold().isin({"nan", "none"}).any():
            continue
        if column != "gene_name":
            adata.var["gene_name"] = labels.to_numpy()
        return "gene_name", column
    labels = pd.Index(adata.var_names.astype(str)).str.strip()
    if len(labels) == 0 or np.any(labels == "") or labels.str.casefold().isin({"nan", "none"}).any():
        raise ValueError("Input has no usable gene_name, symbol column, or var names")
    adata.var["gene_name"] = labels.to_numpy()
    return "gene_name", "var_names"


def _prepare_working_input(raw_path: Path, prepared_path: Path) -> dict[str, Any]:
    """Create a measurement-preserving copy acceptable to the standard profile."""
    import anndata as ad

    raw = ad.read_h5ad(raw_path)
    if raw.n_obs < 3 or raw.n_vars < 3:
        raise ValueError("PBMC3k input must contain at least three cells and genes")
    if not raw.obs_names.is_unique or not raw.var_names.is_unique:
        raise ValueError("Raw H5AD barcodes and feature IDs must be unique for lineage-safe alignment")
    original_shape = tuple(raw.shape)
    original_obs_names = raw.obs_names.astype(str).copy()
    original_var_names = raw.var_names.astype(str).copy()

    if "counts" in raw.layers:
        counts_source = "layers.counts"
        counts = raw.layers["counts"]
    else:
        counts_source = "X"
        counts = raw.X
    _validate_integer_counts(counts, counts_source)
    if tuple(counts.shape) != original_shape:
        raise ValueError(f"Raw count matrix shape {counts.shape} differs from input {original_shape}")
    if "counts" not in raw.layers:
        raw.layers["counts"] = counts.copy() if hasattr(counts, "copy") else np.array(counts, copy=True)

    gene_key, gene_source = _gene_name_column(raw)
    genes = raw.var[gene_key].astype(str)
    matched_markers = sorted(set(genes).intersection({gene for group in PBMC_MARKERS.values() for gene in group}))
    if not matched_markers:
        raise ValueError(
            "None of the fixed PBMC marker genes match the source gene_name values; "
            "supply a raw H5AD with symbols in var['gene_name'] or var names"
        )

    if tuple(raw.shape) != original_shape or not raw.obs_names.astype(str).equals(original_obs_names) or not raw.var_names.astype(str).equals(original_var_names):
        raise AssertionError("Working preparation changed the raw cells or features")
    prepared_path.parent.mkdir(parents=True, exist_ok=True)
    raw.write_h5ad(prepared_path)
    return {
        "n_cells": int(raw.n_obs),
        "n_genes": int(raw.n_vars),
        "shape": list(original_shape),
        "counts_source_before_staging": counts_source,
        "working_counts_layer": "counts",
        "gene_name_column": gene_key,
        "gene_name_source": gene_source,
        "n_unique_gene_names": int(genes.nunique()),
        "matched_fixed_marker_genes": matched_markers,
        "n_obs_columns_before_staging": int(len(raw.obs.columns)),
        "obs_columns_before_staging": [str(value) for value in raw.obs.columns],
        "raw_cells_and_features_preserved_in_working_copy": True,
    }


def _evaluate_reference(clustered: Any, reference: Any, label_key: str = REFERENCE_LABEL_KEY) -> dict[str, Any]:
    """Compare labels only for exact shared barcodes after the EACBP run."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    if "cluster" not in clustered.obs:
        raise ValueError("EACBP result has no obs['cluster'] values")
    if label_key not in reference.obs:
        raise ValueError(f"Reference H5AD has no obs[{label_key!r}] label column")
    if not clustered.obs_names.is_unique or not reference.obs_names.is_unique:
        raise ValueError("Clustered and reference H5AD barcodes must be unique")

    predicted = pd.Series(clustered.obs["cluster"].to_numpy(), index=clustered.obs_names.astype(str))
    truth = pd.Series(reference.obs[label_key].to_numpy(), index=reference.obs_names.astype(str))
    shared = predicted.index.intersection(truth.index, sort=False)
    aligned = pd.DataFrame({"cluster": predicted.loc[shared], "reference": truth.loc[shared]})
    valid = aligned["cluster"].notna() & aligned["reference"].notna()
    aligned = aligned.loc[valid]
    if aligned.empty:
        raise ValueError("No labeled cells have exact barcodes in common with the processed reference")
    aligned["cluster"] = aligned["cluster"].astype(str)
    aligned["reference"] = aligned["reference"].astype(str)
    table = pd.crosstab(aligned["reference"], aligned["cluster"])
    return {
        "status": "computed",
        "label_key": label_key,
        "reference_used_for_pipeline_or_training": False,
        "alignment": "exact obs_names/barcode intersection after clustering",
        "n_reference_cells": int(reference.n_obs),
        "n_clustered_cells": int(clustered.n_obs),
        "n_shared_labeled_cells": int(len(aligned)),
        "n_reference_cells_without_cluster_match": int(reference.n_obs - len(shared)),
        "n_clustered_cells_without_reference_match": int(clustered.n_obs - len(shared)),
        "reference_coverage": {
            "status": "passed" if len(aligned) == reference.n_obs else "failed",
            "criterion": "Every reference barcode has a nonmissing reference label and predicted cluster.",
            "n_reference_cells": int(reference.n_obs),
            "n_shared_barcodes": int(len(shared)),
            "n_shared_labeled_cells": int(len(aligned)),
        },
        "adjusted_rand_index": float(adjusted_rand_score(aligned["reference"], aligned["cluster"])),
        "normalized_mutual_information": float(normalized_mutual_info_score(aligned["reference"], aligned["cluster"])),
        "crosstab": {
            "reference_labels": [str(value) for value in table.index],
            "cluster_labels": [str(value) for value in table.columns],
            "counts": table.to_numpy(dtype=int).tolist(),
        },
    }


def _label_text(value: Any) -> str | None:
    if pd.isna(value):
        return None
    label = str(value).strip()
    if not label or label.casefold() in {"nan", "none", "<na>", "nat"}:
        return None
    return label


def _label_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.casefold())


def _canonical_pbmc_labels() -> dict[str, str]:
    aliases: dict[str, str] = {}
    for canonical, values in PBMC_CELL_TYPE_LABEL_MAP.items():
        for value in values:
            key = _label_key(value)
            if key in aliases and aliases[key] != canonical:
                raise ValueError(f"PBMC label alias {value!r} maps to multiple canonical labels")
            aliases[key] = canonical
    return aliases


def _evaluate_cell_type_annotations(
    clustered: Any,
    reference: Any,
    label_key: str = REFERENCE_LABEL_KEY,
) -> dict[str, Any]:
    """Score marker-derived cell_type labels against fixed human PBMC aliases."""
    from sklearn.metrics import accuracy_score, adjusted_rand_score, f1_score

    if "cell_type" not in clustered.obs:
        raise ValueError("EACBP result has no obs['cell_type'] annotations")
    if label_key not in reference.obs:
        raise ValueError(f"Reference H5AD has no obs[{label_key!r}] label column")
    if not clustered.obs_names.is_unique or not reference.obs_names.is_unique:
        raise ValueError("Clustered and reference H5AD barcodes must be unique")

    predicted = pd.Series(clustered.obs["cell_type"].to_numpy(), index=clustered.obs_names.astype(str))
    truth = pd.Series(reference.obs[label_key].to_numpy(), index=reference.obs_names.astype(str))
    shared = predicted.index.intersection(truth.index, sort=False)
    aligned = pd.DataFrame({"prediction": predicted.loc[shared], "reference": truth.loc[shared]})
    alias_map = _canonical_pbmc_labels()

    raw_prediction = aligned["prediction"].map(_label_text)
    raw_reference = aligned["reference"].map(_label_text)
    canonical_prediction = raw_prediction.map(
        lambda value: alias_map.get(_label_key(value)) if value is not None else None
    )
    canonical_reference = raw_reference.map(
        lambda value: alias_map.get(_label_key(value)) if value is not None else None
    )

    reference_missing = raw_reference.isna()
    reference_unmapped = raw_reference.notna() & canonical_reference.isna()
    prediction_missing = raw_prediction.isna()
    prediction_unmapped = raw_prediction.notna() & canonical_prediction.isna()
    scored = canonical_reference.notna()

    def display_labels(raw_labels, canonical_labels, missing_mask, unmapped_mask, role):
        values = []
        for raw, canonical, is_missing, is_unmapped in zip(
            raw_labels, canonical_labels, missing_mask, unmapped_mask
        ):
            if is_missing:
                values.append(f"<missing {role}>")
            elif is_unmapped:
                values.append(f"<unmapped {role}: {raw}>")
            else:
                values.append(str(canonical))
        return values

    table = pd.crosstab(
        pd.Series(display_labels(raw_reference, canonical_reference, reference_missing, reference_unmapped, "reference"), index=shared),
        pd.Series(display_labels(raw_prediction, canonical_prediction, prediction_missing, prediction_unmapped, "prediction"), index=shared),
    )
    common = {
        "label_key": label_key,
        "label_map": {canonical: list(aliases) for canonical, aliases in PBMC_CELL_TYPE_LABEL_MAP.items()},
        "reference_used_for_pipeline_or_training": False,
        "alignment": "exact obs_names/barcode intersection after clustering",
        "n_clustered_cells": int(clustered.n_obs),
        "n_reference_cells": int(reference.n_obs),
        "n_shared_cells": int(len(shared)),
        "n_reference_labels_missing": int(reference_missing.sum()),
        "n_reference_labels_unmapped": int(reference_unmapped.sum()),
        "n_prediction_labels_missing": int(prediction_missing.sum()),
        "n_prediction_labels_unmapped": int(prediction_unmapped.sum()),
        "reference_unmapped_label_counts": raw_reference.loc[reference_unmapped].value_counts().to_dict(),
        "prediction_unmapped_label_counts": raw_prediction.loc[prediction_unmapped].value_counts().to_dict(),
        "annotation_crosstab": {
            "reference_labels": [str(value) for value in table.index],
            "prediction_labels": [str(value) for value in table.columns],
            "counts": table.to_numpy(dtype=int).tolist(),
        },
    }
    if not scored.any():
        return {
            **common,
            "status": "not_estimable",
            "metrics_scope": "No exact-barcode reference labels matched the fixed PBMC alias map.",
            "n_cells_scored": 0,
            "scored_reference_labels": [],
            "accuracy": None,
            "macro_f1": None,
            "annotation_ari": None,
        }

    true_labels = canonical_reference.loc[scored].astype(str).to_numpy()
    prediction_labels = canonical_prediction.loc[scored].astype(object).to_numpy()
    prediction_values = [
        str(value) if value is not None and not pd.isna(value)
        else f"<unmapped prediction: {raw_prediction.loc[index]}>"
        if pd.notna(raw_prediction.loc[index])
        else "<missing prediction>"
        for index, value in zip(canonical_prediction.loc[scored].index, prediction_labels)
    ]
    observed_reference_classes = sorted(set(true_labels))
    return {
        **common,
        "status": "computed",
        "metrics_scope": "Exact shared barcodes with a nonmissing reference label mapped by the fixed alias map; unknown predictions count as errors.",
        "n_cells_scored": int(scored.sum()),
        "scored_reference_labels": observed_reference_classes,
        "accuracy": float(accuracy_score(true_labels, prediction_values)),
        "macro_f1": float(f1_score(
            true_labels, prediction_values, labels=observed_reference_classes,
            average="macro", zero_division=0,
        )),
        "annotation_ari": float(adjusted_rand_score(true_labels, prediction_values)),
    }


def evaluate_saved_annotations(clustered_path: Path, reference_path: Path) -> dict[str, Any]:
    """Evaluate an already-saved cluster H5AD without rerunning the pipeline."""
    import anndata as ad

    clustered = ad.read_h5ad(clustered_path, backed="r")
    reference = ad.read_h5ad(reference_path, backed="r")
    try:
        return _evaluate_cell_type_annotations(clustered, reference)
    finally:
        clustered.file.close()
        reference.file.close()


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(_jsonable(payload), ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _validate_external_work_dir_path(work_dir: PurePath, project_root: PurePath) -> PurePath:
    """Apply portable path-policy checks without touching the filesystem."""
    if not work_dir.is_absolute():
        raise ValueError("--work-dir must be an absolute path on the external validation volume")
    try:
        work_dir.relative_to(project_root)
    except ValueError:
        pass
    else:
        raise ValueError("--work-dir must be outside the repository/workspace so generated data stays on the validation volume")
    if isinstance(work_dir, PureWindowsPath) and work_dir.drive.casefold() == "c:":
        raise ValueError("--work-dir must not use the C: drive; choose an external validation volume such as F:")
    return work_dir


def _validate_external_work_dir(work_dir: Path, *, project_root: Path = PROJECT_ROOT) -> Path:
    _validate_external_work_dir_path(work_dir, project_root)
    resolved = work_dir.expanduser().resolve()
    _validate_external_work_dir_path(resolved, project_root.resolve())
    return resolved


def run_validation(raw_path: Path, reference_path: Path, work_dir: Path) -> tuple[dict[str, Any], Path]:
    """Run the standard profile and write a JSON receipt plus clustered H5AD."""
    import anndata as ad

    raw_path = raw_path.expanduser().resolve()
    reference_path = reference_path.expanduser().resolve()
    work_dir = _validate_external_work_dir(work_dir)
    if not raw_path.is_file() or not reference_path.is_file():
        raise FileNotFoundError("Both --data and --reference must be existing H5AD files")
    if raw_path == reference_path:
        raise ValueError("--data and --reference must be different files")
    work_dir.mkdir(parents=True, exist_ok=True)
    run_root = work_dir / f"pbmc3k-validation-{uuid4().hex}"
    run_root.mkdir()
    raw_hash_before = _sha256(raw_path)
    reference_hash_before = _sha256(reference_path)
    prepared_path = run_root / "input" / "pbmc3k_raw_counts.h5ad"
    result_json = run_root / "public_pbmc_validation.json"
    result_h5ad = run_root / "pbmc3k_standard_clusters.h5ad"
    result: dict[str, Any] = {
        "schema_version": 1,
        "dataset": "PBMC3k",
        "analysis": "standard_profile_unsupervised_clustering",
        "status": "failed",
        "statuses": {"engineering": "not_run", "benchmark": "not_run"},
        "benchmark_scope": "Input integrity and complete reference comparison; metrics are descriptive.",
        "scientific_interpretation": "pending_domain_review",
        "biological_validity_claimed": False,
        "source": {"path": str(raw_path), "sha256_before": raw_hash_before,
                   "expected_sha256": PUBLIC_INPUT_SHA256["pbmc3k_raw"]},
        "reference": {"path": str(reference_path), "label_key": REFERENCE_LABEL_KEY,
                       "used_for_pipeline_or_training": False,
                       "sha256_before": reference_hash_before,
                       "expected_sha256": PUBLIC_INPUT_SHA256["pbmc3k_processed"]},
        "run_root": str(run_root),
        "outputs": {"result_json": str(result_json), "clustered_h5ad": None},
        "statistical_limits": {
            "differential_expression": "not_applicable",
            "reason": "PBMC3k input has no independent donor/condition design; no donor-level inference is requested or synthesized.",
        },
    }

    try:
        for name in ("source", "reference"):
            record = result[name]
            record["matches_pinned_sha256"] = record["sha256_before"] == record["expected_sha256"]
            if not record["matches_pinned_sha256"]:
                raise ValueError(
                    f"Public PBMC3k {name} SHA256 mismatch: expected {record['expected_sha256']}, "
                    f"got {record['sha256_before']}"
                )
        staged = _prepare_working_input(raw_path, prepared_path)
        result["input"] = {**staged, "working_h5ad": str(prepared_path),
                           "working_h5ad_sha256": _sha256(prepared_path)}
        result["method"] = {
            "method_profile": "standard",
            "integration": "no_correction_v1",
            "clustering": "scanpy_leiden_umap_v1",
            "qc": {"min_genes": 10, "max_mito_pct": 20.0,
                   "mitochondrial_metric_source_is_recorded_in_qc_artifact": True},
            "marker_reference": PBMC_MARKERS,
            "reference_label_key": REFERENCE_LABEL_KEY,
        }

        from eacbp.application.study_service import run_study
        from eacbp.artifact.registry import ArtifactRegistry
        from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest

        study_id = "public_pbmc3k"
        manifest = StudyManifest(
            study_id=study_id,
            biological_design=BiologicalDesign(species="human", tissue="PBMC", target_cell_types=[]),
            data=DataSpec(),
        )
        config = {
            "mode": "real",
            "method_profile": "standard",
            "method_overrides": {"integration": "no_correction_v1"},
            "capability_parameters": {
                "qc": {"min_genes": 10, "max_mito_pct": 20.0},
                "clustering": {"marker_reference": PBMC_MARKERS},
            },
        }
        run_dir = run_root / "eacbp-run"
        result["statuses"]["engineering"] = "running"
        summary = run_study(manifest=manifest, config=config, data=prepared_path, run_dir=run_dir)
        result["run"] = _jsonable(summary)
        result["statuses"]["engineering"] = "passed" if summary.get("status") == "success" else "failed"

        registry = ArtifactRegistry(str(run_dir / "artifacts"))
        qc_uri = f"adata://{study_id}/qc/v1"
        if registry.exists(qc_uri):
            result["qc"] = _jsonable(registry.get_metadata(qc_uri).summary_metrics)
        artifact_uri = f"adata://{study_id}/annotated/v4"
        result["cluster_artifact_uri"] = artifact_uri
        if registry.exists(artifact_uri):
            _, clustered_data = registry.get_scdata(artifact_uri)
            clustered = clustered_data.to_anndata()
            clustered.write_h5ad(result_h5ad)
            result["outputs"]["clustered_h5ad"] = str(result_h5ad)
            result["cluster_result"] = {
                "n_cells": int(clustered.n_obs),
                "n_genes": int(clustered.n_vars),
                "n_clusters": int(clustered.obs["cluster"].nunique()) if "cluster" in clustered.obs else None,
                "cell_type_counts": clustered.obs["cell_type"].astype(str).value_counts().to_dict()
                if "cell_type" in clustered.obs else {},
            }
            reference = ad.read_h5ad(reference_path, backed="r")
            try:
                result["reference_evaluation"] = _evaluate_reference(clustered, reference)
                annotations = _evaluate_cell_type_annotations(clustered, reference)
                result["reference_evaluation"]["annotation_evaluation"] = annotations
                annotation_complete = (
                    annotations["status"] == "computed"
                    and annotations["n_cells_scored"] == reference.n_obs
                )
                result["reference_evaluation"]["annotation_coverage"] = {
                    "status": "passed" if annotation_complete else "failed",
                    "criterion": "Every reference barcode has a reference label scored by the fixed PBMC alias map.",
                    "n_reference_cells": int(reference.n_obs),
                    "n_cells_scored": annotations["n_cells_scored"],
                }
                comparison_complete = (
                    result["reference_evaluation"]["reference_coverage"]["status"] == "passed"
                    and annotation_complete
                )
                result["statuses"]["benchmark"] = "passed" if comparison_complete else "failed"
                if not comparison_complete:
                    result["failure_note"] = "Reference comparison did not cover every reference cell and fixed reference label."
            finally:
                reference.file.close()
        else:
            result["statuses"]["benchmark"] = "failed"
            result["reference_evaluation"] = {
                "status": "skipped",
                "reason": "EACBP did not register the audited clustering artifact; reference labels were not read.",
            }

        if result["statuses"]["engineering"] != "passed":
            result["failure_note"] = "Run or audit failures are retained as reported; comparison metrics, if present, are diagnostic only."
    except Exception as exc:
        if result["statuses"]["engineering"] == "running":
            result["statuses"]["engineering"] = "failed"
        result["statuses"]["benchmark"] = "failed"
        result["error"] = {"type": type(exc).__name__, "message": str(exc)}
    finally:
        integrity_failures = []
        for name, path in (("source", raw_path), ("reference", reference_path)):
            record = result[name]
            try:
                record["sha256_after"] = _sha256(path)
                record["unchanged"] = record["sha256_before"] == record["sha256_after"]
            except OSError as exc:
                record["sha256_after"] = None
                record["unchanged"] = False
                record["integrity_error"] = f"{type(exc).__name__}: {exc}"
            if not record["unchanged"]:
                integrity_failures.append(name)
        if integrity_failures:
            result["statuses"]["benchmark"] = "failed"
            result["integrity_failures"] = integrity_failures
            result["failure_note"] = "Public input SHA256 changed or became unavailable during the run."
    result["status"] = "success" if all(value == "passed" for value in result["statuses"].values()) else "failed"
    _write_json(result_json, result)
    return result, result_json


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="read-only raw PBMC3k H5AD")
    parser.add_argument("--reference", type=Path, required=True, help="processed reference H5AD; labels are used only after clustering")
    parser.add_argument("--work-dir", type=Path, required=True, help="required output directory, preferably on the validation drive")
    args = parser.parse_args(argv)
    result, result_path = run_validation(args.data, args.reference, args.work_dir)
    print(json.dumps({"status": result["status"], "statuses": result["statuses"], "result_json": str(result_path),
                      "clustered_h5ad": result.get("outputs", {}).get("clustered_h5ad"),
                      "reference_evaluation": result.get("reference_evaluation")},
                     ensure_ascii=False, indent=2, allow_nan=False))
    return 0 if result["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
