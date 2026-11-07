"""Durable, task-scoped receipts for actual PyDESeq2 fit attempts.

Receipts are diagnostic records, not statistical evidence of robustness. The
caller publishes the batch as a registered JSON task output so the ordinary
transaction, hashing, and independent admission also cover these records.
"""
from __future__ import annotations

from contextvars import ContextVar
from datetime import date, datetime, timezone
import hashlib
import importlib.metadata
import json
import platform
from typing import Any, Mapping
import warnings

import numpy as np
import pandas as pd

from eacbp.artifact.uri import ArtifactURI


SCHEMA = "eacbp.pydeseq2-fit-receipts.v1"
_active: ContextVar[FitReceipt | None] = ContextVar("eacbp_fit_receipt", default=None)


def fit_receipts_uri(primary_uri: str) -> str:
    uri = ArtifactURI.parse(primary_uri)
    return ArtifactURI("json", uri.study_id, uri.name + "_fit_receipts", uri.version).to_string()


def _plain(value: Any) -> Any:
    if isinstance(value, np.generic):
        return _plain(value.item())
    if isinstance(value, (datetime, date, pd.Timestamp)):
        return value.isoformat()
    if isinstance(value, pd.Series):
        return _plain(value.to_list())
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(item) for item in value]
    return value


def canonical_sha256(value: Any) -> str:
    encoded = json.dumps(_plain(value), sort_keys=True, ensure_ascii=False,
                         separators=(",", ":"), allow_nan=False).encode("utf-8")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def frame_record(frame: pd.DataFrame) -> dict[str, Any]:
    record = _plain(frame.to_dict(orient="split"))
    # String identities are explicit; numeric count and design entries retain
    # their original values. There is no replacement of statistical NaN/Inf.
    record["index"] = [str(value) for value in frame.index]
    record["columns"] = [str(value) for value in frame.columns]
    return {"frame": record, "sha256": canonical_sha256(record)}


def software_versions() -> dict[str, str]:
    result = {"python": platform.python_version()}
    for package in ("pydeseq2", "numpy", "pandas", "scipy", "formulaic", "anndata"):
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = "unavailable"
    return result


class FitReceipt:
    """Capture one attempt, including an exception or a successful fallback."""

    def __init__(self, *, task_id: str, method: str, counts: pd.DataFrame,
                 metadata: pd.DataFrame, design_formula: str,
                 contrast: Mapping[str, Any], params: Mapping[str, Any],
                 left_out_donor: str | None):
        count_record = frame_record(counts)
        metadata_record = frame_record(metadata)
        identity = {"task_id": task_id, "method": method,
                    "left_out_donor": left_out_donor,
                    "counts_sha256": count_record["sha256"],
                    "metadata_sha256": metadata_record["sha256"],
                    "design_formula": design_formula, "contrast": dict(contrast)}
        self.record: dict[str, Any] = {
            **identity,
            "fit_id": canonical_sha256(identity),
            "fit_role": "leave_one_donor_out" if left_out_donor is not None else "full_model",
            "sample_ids": list(count_record["frame"]["index"]),
            "counts": count_record, "metadata": metadata_record,
            "requested_options": {key: _plain(params[key]) for key in (
                "fit_type", "size_factors_fit_type", "min_replicates", "low_memory",
                "n_cpus", "n_processes", "alpha", "cooks_filter", "independent_filter"
            ) if key in params},
            "effective_options": {}, "design_matrix": None,
            "actual_dispersion_trend": None, "dispersion_trend_source": None,
            "software_versions": software_versions(),
            "warnings": [],
            "warning_capture_scope": "Python warnings emitted in the fitting process; worker-process stderr is not captured",
            "status": "attempted", "error": None,
            "started_at": datetime.now(timezone.utc).isoformat(),
        }
        self.dds = None

    def __enter__(self):
        self._token = _active.set(self)
        self._warning_context = warnings.catch_warnings(record=True)
        self._captured = self._warning_context.__enter__()
        # Capture warnings per attempt, including repeats from another fold.
        # Existing filters still apply: recording must not turn a warning that
        # the caller treats as an error into a successful fit.
        return self

    def observe_dds(self, dds: Any, effective_options: Mapping[str, Any]) -> None:
        self.dds = dds
        self.record["effective_options"].update(_plain(effective_options))
        design = dds.obsm.get("design_matrix")
        if isinstance(design, pd.DataFrame):
            self.record["design_matrix"] = frame_record(design)
        observed = {key: getattr(dds, key) for key in (
            "fit_type", "size_factors_fit_type", "min_replicates", "low_memory",
            "min_mu", "min_disp", "max_disp", "refit_cooks", "beta_tol"
        ) if hasattr(dds, key)}
        self.record["effective_options"].update(_plain(observed))

    def observe_stats(self, stats_obj: Any) -> None:
        self.record["effective_options"].update({key: _plain(getattr(stats_obj, key))
            for key in ("alpha", "cooks_filter", "independent_filter", "lfc_null", "alt_hypothesis")
            if hasattr(stats_obj, key)})
        vector = getattr(stats_obj, "contrast_vector", None)
        if vector is not None:
            self.record["contrast_vector"] = _plain(vector)

    def __exit__(self, exc_type, exc, traceback):
        try:
            self.record["warnings"] = [{"message": str(item.message),
                "category": item.category.__name__, "filename": str(item.filename),
                "lineno": int(item.lineno)} for item in self._captured]
            if self.dds is not None:
                # fit_type is the requested option. Only this fitted-state
                # field establishes which trend the installed library used.
                actual = getattr(self.dds, "uns", {}).get("disp_function_type")
                if actual in {"parametric", "mean"}:
                    self.record["actual_dispersion_trend"] = actual
                    self.record["dispersion_trend_source"] = "dds.uns.disp_function_type"
            self.record["status"] = "failed" if exc is not None else "completed"
            self.record["error"] = None if exc is None else {"type": exc_type.__name__, "message": str(exc)}
            self.record["finished_at"] = datetime.now(timezone.utc).isoformat()
            self.record["receipt_sha256"] = canonical_sha256(self.record)
        finally:
            self._warning_context.__exit__(exc_type, exc, traceback)
            _active.reset(self._token)
        return False


def observe_dds(dds: Any, effective_options: Mapping[str, Any]) -> None:
    receipt = _active.get()
    if receipt is not None:
        receipt.observe_dds(dds, effective_options)


def observe_stats(stats_obj: Any) -> None:
    receipt = _active.get()
    if receipt is not None:
        receipt.observe_stats(stats_obj)
