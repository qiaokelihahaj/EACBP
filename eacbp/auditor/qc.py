"""Independently verify native QC dispositions against the parent cell matrix."""

import numpy as np
import pandas as pd

from eacbp.auditor.base import ValidationCheck, ValidationSeverity


def validate_qc_disposition(contract, result, registry, output_meta, output):
    """Allow undeclared legacy records; declared new records must be complete."""
    qc = output.uns.get("qc", {})
    if not isinstance(qc, dict):
        return ValidationCheck(check_name="qc_cell_disposition_integrity", passed=False,
                               severity=ValidationSeverity.ERROR,
                               message="QC provenance must be a structured record.",
                               metrics={"cell_disposition_available": False, "all_input_cells_assessed": False})
    declaration = [qc, output_meta.summary_metrics, result.metrics]
    declared = any("cell_disposition_schema_version" in item for item in declaration)
    present = "cell_disposition" in qc
    if not declared and not present:
        return ValidationCheck(
            check_name="qc_cell_disposition_integrity", passed=True,
            severity=ValidationSeverity.INFO,
            message="Legacy QC artifact has no cell disposition table; per-cell disposition was not assessed.",
            metrics={"cell_disposition_available": False, "all_input_cells_assessed": False},
        )
    metrics = {"cell_disposition_available": present, "all_input_cells_assessed": False}
    try:
        def require(condition, message):
            if not condition:
                raise ValueError(message)

        require(present, "Declared native QC disposition table is missing.")
        require(all(item.get("cell_disposition_schema_version") == 1 for item in declaration),
                "QC disposition schema declarations disagree or are incomplete.")
        require(len(contract.input_artifacts) == 1, "QC must have one parent input.")
        input_uri = contract.input_artifacts[0]
        input_meta, source = registry.get_scdata(input_uri)
        require(output_meta.parent_uris == [input_uri], "QC lineage does not match the contracted input.")
        require(qc.get("input_artifact_uri") == input_uri and qc.get("input_sha256") == input_meta.sha256_hash,
                "QC disposition table is not bound to the verified input hash.")
        ids = source.obs.index.astype(str).to_numpy()
        require(pd.Index(ids).is_unique, "Input cell identifiers are not unique.")
        counts = source.layers.get("counts", source.X)
        values = counts.data if hasattr(counts, "tocsr") else np.asarray(counts)
        require(np.isfinite(values).all() and (values >= 0).all(), "Input counts are not finite and nonnegative.")
        genes = np.asarray((counts > 0).sum(axis=1)).ravel()
        totals = np.asarray(counts.sum(axis=1)).ravel()
        minimum = float(contract.parameters.get("min_genes", 100))
        maximum = float(contract.parameters.get("max_mito_pct", 20.0))
        require(np.isfinite(minimum) and minimum >= 0 and np.isfinite(maximum) and 0 <= maximum <= 100,
                "QC thresholds are invalid.")
        require(qc.get("min_genes") == minimum and qc.get("max_mito_pct") == maximum,
                "Recorded QC thresholds differ from the contract.")
        require(output_meta.parameters == {"min_genes": minimum, "max_mito_pct": maximum},
                "QC artifact parameters differ from the contract.")
        mito = None
        metric_source = "unavailable"
        n_mito = None
        if "percent_mito" in source.obs:
            mito = pd.to_numeric(source.obs["percent_mito"], errors="raise").to_numpy(dtype=float)
            metric_source = "obs.percent_mito"
        else:
            if "mt" in source.var:
                require(source.var["mt"].dtype == bool and not source.var["mt"].isna().any(),
                        "Input mitochondrial annotations are not complete booleans.")
                mt = source.var["mt"].to_numpy(dtype=bool)
                candidate_source = "var.mt"
            else:
                symbols = source.var["gene_name"] if "gene_name" in source.var else source.var.index
                mt = np.asarray(pd.Index(symbols).astype(str).str.lower().str.startswith("mt-"))
                candidate_source = "gene_symbol_prefix"
            n_mito = int(mt.sum())
            if n_mito:
                mitochondrial_counts = np.asarray(counts[:, mt].sum(axis=1)).ravel()
                mito = np.divide(mitochondrial_counts * 100.0, totals,
                                 out=np.zeros(source.n_obs, dtype=float), where=totals > 0)
                metric_source = candidate_source
        if mito is not None:
            require(np.isfinite(mito).all() and ((mito >= 0) & (mito <= 100)).all(),
                    "Input mitochondrial percentages are invalid.")
        low = genes < minimum
        high = mito > maximum if mito is not None else np.zeros(source.n_obs, dtype=bool)
        retained = ~(low | high)
        available = np.full(source.n_obs, mito is not None, dtype=bool)
        if metric_source in ("var.mt", "gene_symbol_prefix"):
            available &= totals > 0
        mito_record = np.full(source.n_obs, np.nan, dtype=float)
        if mito is not None:
            mito_record[available] = mito[available]
        reasons = np.full(source.n_obs, "retained", dtype=object)
        reasons[low] = "low_gene_count"
        reasons[high] = "high_mito_percentage"
        reasons[low & high] = "low_gene_count_and_high_mito_percentage"
        expected = {"ids": ids, "n_genes": genes, "total_counts": totals, "mito_pct": mito_record,
                    "mito_pct_available": available, "low_gene_flag": low, "high_mito_flag": high,
                    "retained": retained, "reason": reasons}
        table = qc["cell_disposition"]
        require(isinstance(table, dict) and set(table) == set(expected), "QC disposition columns are incomplete or unexpected.")
        for key, expected_values in expected.items():
            actual = np.asarray(table[key])
            require(actual.shape == (source.n_obs,), f"QC disposition column {key} does not cover every input cell.")
            if expected_values.dtype == bool:
                require(actual.dtype == bool, f"QC disposition column {key} is not boolean.")
            equal = np.array_equal(actual, expected_values, equal_nan=True) if key == "mito_pct" else np.array_equal(actual, expected_values)
            require(equal, f"QC disposition column {key} disagrees with independently recomputed input values.")
        expected_counts = {name: int(np.count_nonzero(reasons == name)) for name in
                           ("retained", "low_gene_count", "high_mito_percentage", "low_gene_count_and_high_mito_percentage")}
        for record in declaration:
            require(record.get("cell_disposition_counts") == expected_counts, "QC disposition counts disagree with input cells.")
            warning = record.get("interpretation_warning", "")
            require("does not independently establish cellular damage" in warning and "not a damage label" in warning,
                    "QC disposition interpretation warning is missing.")
        require(qc.get("cell_disposition_id_source") == "obs.index", "QC disposition identifier source is invalid.")
        require(qc.get("counts_source") == ("layers.counts" if "counts" in source.layers else "X") and
                qc.get("mitochondrial_metric_source") == metric_source and
                qc.get("mitochondrial_filtering_applied") == (mito is not None) and
                qc.get("n_mitochondrial_genes") == n_mito, "QC metric provenance disagrees with the input.")
        require(np.array_equal(output.obs.index.astype(str).to_numpy(), ids[retained]),
                "QC output cell identifiers/order disagree with the threshold mask.")
        require(output.var.index.equals(source.var.index), "QC unexpectedly changed the gene axis.")
        def same_matrix(actual, expected_matrix):
            if hasattr(actual, "tocsr") and hasattr(expected_matrix, "tocsr"):
                return actual.shape == expected_matrix.shape and (actual != expected_matrix).nnz == 0
            if hasattr(actual, "tocsr") or hasattr(expected_matrix, "tocsr"):
                return False
            return np.array_equal(actual, expected_matrix)
        require(same_matrix(output.X, source.X[retained]), "QC output expression values disagree with the retained input.")
        require(same_matrix(output.layers.get("counts", output.X), counts[retained]),
                "QC output counts disagree with the retained input.")
        for key, expected_values in (("qc_passed", np.ones(int(retained.sum()), dtype=bool)),
                                     ("n_genes_by_counts", genes[retained]), ("total_counts", totals[retained])):
            require(key in output.obs and np.array_equal(output.obs[key].to_numpy(), expected_values),
                    f"QC retained-cell annotation {key} disagrees with the input.")
        if mito is not None:
            require("percent_mito" in output.obs and np.allclose(output.obs["percent_mito"], mito[retained]),
                    "QC retained mitochondrial percentages disagree with the input.")
        metrics.update({"all_input_cells_assessed": True, "initial_cells": source.n_obs,
                        "retained_cells": int(retained.sum()), "cell_disposition_counts": expected_counts})
        return ValidationCheck(check_name="qc_cell_disposition_integrity", passed=True,
                               severity=ValidationSeverity.INFO,
                               message="Independently recomputed all input-cell QC metrics, dispositions and the retained matrix subset.",
                               metrics=metrics)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as exc:
        return ValidationCheck(check_name="qc_cell_disposition_integrity", passed=False,
                               severity=ValidationSeverity.ERROR, message=str(exc), metrics=metrics)
