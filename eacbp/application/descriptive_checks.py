"""Recompute selected raw-UMI summaries and check external descriptive tables.

No cell-level hypothesis tests or module rescoring occur here. Cell identities
bind mask comparisons; matching cell counts alone never establish equivalence.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

import numpy as np
import pandas as pd


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _columns(frame, names, label):
    _require(frame.columns.is_unique, f"{label}: duplicate column names")
    missing = sorted(set(names) - set(frame.columns))
    _require(not missing, f"{label}: missing required columns {missing}")


def _labels(series, label):
    _require(not series.isna().any(), f"{label}: missing identity")
    result = series.astype(str)
    _require(bool(result.str.strip().eq(result).all() and result.ne("").all()),
             f"{label}: empty or padded identity")
    _require(not result.str.casefold().isin({"nan", "none", "null", "<na>"}).any(),
             f"{label}: invalid identity")
    return result


def _numbers(series, label, *, integer=False, nonnegative=False, nullable=False):
    _require(not series.map(lambda value: isinstance(value, (bool, np.bool_))).any(),
             f"{label}: booleans are not numeric counts")
    try:
        values = pd.to_numeric(series, errors="raise").to_numpy(dtype=float, na_value=np.nan)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{label}: expected numeric values") from exc
    valid = np.isfinite(values)
    _require(bool((valid | np.isnan(values)).all()) if nullable else bool(valid.all()),
             f"{label}: values must be finite{' or missing' if nullable else ''}")
    if nonnegative:
        _require(bool((values[valid] >= 0).all()), f"{label}: values must be nonnegative")
    if integer:
        _require(bool((values[valid] == np.floor(values[valid])).all()),
                 f"{label}: values must be integers")
        _require(bool((values[valid] <= 2**53 - 1).all()), f"{label}: integer exceeds exact numeric range")
    return values


def _booleans(series, label):
    def parse(value):
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if isinstance(value, str) and value in {"True", "False", "true", "false"}:
            return value.casefold() == "true"
        raise ValueError(f"{label}: expected strict boolean True/False, found {value!r}")
    return series.map(parse).to_numpy(dtype=bool)


def _equal(actual, expected, label):
    if expected is None:
        _require(pd.isna(actual), f"{label}: expected missing value")
        return
    if isinstance(expected, str):
        _require(actual == expected, f"{label}: expected {expected!r}, found {actual!r}")
        return
    value = _numbers(pd.Series([actual]), label)[0]
    if isinstance(expected, (int, np.integer)):
        _require(value == expected, f"{label}: expected exact integer {expected!r}, found {actual!r}")
        return
    _require(bool(np.isclose(value, expected, rtol=1e-9, atol=1e-10)),
             f"{label}: expected {expected!r}, found {actual!r}")


def _digest(ids):
    return hashlib.sha256(json.dumps(sorted(ids), ensure_ascii=False,
                                    separators=(",", ":")).encode()).hexdigest()


def _sign_status(values):
    if len(values) < 2:
        return "insufficient_distinct_comparable_cell_sets"
    signs = set(np.sign(values))
    if signs == {1}:
        return "consistent_positive"
    if signs == {-1}:
        return "consistent_negative"
    if signs == {0}:
        return "consistent_zero"
    return "sign_flip" if {-1, 1} <= signs else "mixed_zero_and_nonzero"


def check_descriptive_tables(cells: pd.DataFrame, modules: pd.DataFrame, genes: pd.DataFrame, *,
                             reference_sample: str, query_sample: str,
                             population: str, min_cells: int) -> dict[str, Any]:
    """Check project-format tables; raise on missing, contradictory or unsafe data.

    ``population`` is the reviewed subtype label covering the supplied cell
    universe, when a subtype summary is present. Numerical clusters are treated
    only as existing expression states. The full-feature per-cell total_counts
    is a supplied denominator, not reconstructed from the selected gene panel.
    """
    _require(isinstance(min_cells, int) and not isinstance(min_cells, bool) and min_cells >= 1,
             "min_cells must be a positive integer descriptive threshold")
    _labels(pd.Series([reference_sample, query_sample, population]), "requested labels")
    _require(reference_sample != query_sample, "reference_sample and query_sample must differ")
    cells, modules, genes = cells.copy(), modules.copy(), genes.copy()
    _columns(cells, ["cell_id", "sample_id", "cluster", "total_counts"], "cells")
    _require(not cells.empty, "cells: input is empty")
    for column in ("cell_id", "sample_id", "cluster"):
        cells[column] = _labels(cells[column], f"cells.{column}")
    _require(not cells.cell_id.duplicated().any(), "cells: duplicate cell_id")
    _require(set(cells.sample_id) == {reference_sample, query_sample},
             "cells: sample_id must contain exactly the reference and query libraries")
    totals = _numbers(cells.total_counts, "cells.total_counts", integer=True, nonnegative=True)
    _require(sum(int(value) for value in totals) <= 2**53 - 1,
             "cells: aggregate total_counts exceeds exact numeric range")
    raw_columns = [name for name in cells if name.endswith("_raw_UMI")]
    _require(bool(raw_columns), "cells: no selected gene *_raw_UMI columns")
    raw_genes = [name[:-8] for name in raw_columns]
    _labels(pd.Series(raw_genes), "selected genes")
    raw = np.column_stack([_numbers(cells[name], f"cells.{name}", integer=True, nonnegative=True)
                           for name in raw_columns])
    _require(bool((raw.sum(axis=1) <= totals).all()),
             "cells: selected gene raw UMI sum exceeds full-library total_counts")
    mask_names = [name[5:-5] for name in cells if name.startswith("mask_") and name.endswith("_pass")]
    known_names = [name[5:-15] for name in cells if name.startswith("mask_") and name.endswith("_metadata_known")]
    _require(set(mask_names) == set(known_names) and "all" in mask_names,
             "cells: every mask requires pass and metadata_known columns, including all")
    mask_names = ["all", *[name for name in mask_names if name != "all"]]
    masks = {}
    for name in mask_names:
        passed = _booleans(cells[f"mask_{name}_pass"], f"mask_{name}_pass")
        known = _booleans(cells[f"mask_{name}_metadata_known"], f"mask_{name}_metadata_known")
        _require(not (passed & ~known).any(), f"mask {name}: missing metadata cannot pass")
        masks[name] = (passed, known)
    _require(bool(masks["all"][0].all() and masks["all"][1].all()), "mask all must include every cell with known membership")
    all_cells = np.ones(len(cells), dtype=bool)
    groups = {("pooled", "ALL_REVIEWED_EPITHELIAL"): all_cells}
    for state in sorted(set(cells.cluster)):
        groups[("subcluster_state", "subcluster_" + state)] = cells.cluster.eq(state).to_numpy()
    gene_groups = dict(groups)
    module_keys = ["group_scope", "group_label", "mask", "module_id"]
    side_fields = ["n_cells", "n_finite_scores", "median_score", "mask_status", "n_metadata_missing"]
    _columns(modules, [*module_keys, "reference_sample", "query_sample", "comparison_status", "median_score_delta",
                       *[f"{field}_{side}" for field in side_fields for side in ("reference", "query")]], "modules")
    _require(not modules.empty, "modules: input is empty")
    for name in [*module_keys, "reference_sample", "query_sample", "comparison_status"]:
        modules[name] = _labels(modules[name], f"modules.{name}")
    _require(not modules.duplicated(module_keys).any(), "modules: duplicate group/mask/module row")
    _require(modules.reference_sample.eq(reference_sample).all() and modules.query_sample.eq(query_sample).all(),
             "modules: reference/query orientation disagrees with requested libraries")
    if modules.group_scope.eq("subtype").any():
        labels = set(modules.loc[modules.group_scope.eq("subtype"), "group_label"])
        _require(labels == {population}, "modules: subtype label must equal the declared whole-cell population")
        groups[("subtype", population)] = all_cells
        if "group_scope" in genes and genes.group_scope.eq("subtype").any():
            gene_groups[("subtype", population)] = all_cells
    _require(set(map(tuple, modules[["group_scope", "group_label"]].to_numpy())) == set(groups),
             "modules: required pooled/state groups missing or unknown group labels present")
    module_ids = sorted(set(modules.module_id))
    expected_module_keys = {(scope, label, mask, module) for scope, label in groups
                            for mask in mask_names for module in module_ids}
    _require(set(map(tuple, modules[module_keys].to_numpy())) == expected_module_keys,
             "modules: missing or unexpected required group/mask/module rows")
    for side in ("reference", "query"):
        for field in ("n_cells", "n_finite_scores", "n_metadata_missing"):
            _numbers(modules[f"{field}_{side}"], f"modules.{field}_{side}", integer=True, nonnegative=True)
        _numbers(modules[f"median_score_{side}"], f"modules.median_score_{side}", nullable=True)
    _numbers(modules.median_score_delta, "modules.median_score_delta", nullable=True)

    selections = {}
    for group, members in groups.items():
        for mask, (passed, known) in masks.items():
            for sample in (reference_sample, query_sample):
                base = members & cells.sample_id.eq(sample).to_numpy()
                n_input, n_known = int(base.sum()), int((base & known).sum())
                available = n_input == 0 or n_known > 0
                selected = base & passed
                status = ("unavailable" if not available else "missing_cells" if n_input == 0 else
                          "available_with_missing_metadata" if n_known < n_input else "available")
                selections[(*group, mask, sample)] = {
                    "selected": selected, "n_input_cells": n_input, "n_metadata_known": n_known,
                    "n_metadata_missing": n_input - n_known, "n_cells": int(selected.sum()),
                    "mask_available": available, "mask_status": status,
                    "total_selected_library_UMI": int(totals[selected].sum()) if available else None}

    gene_keys = ["group_scope", "group_label", "mask", "sample_id", "gene_symbol"]
    gene_metrics = ["n_cells", "total_selected_library_UMI", "raw_count_sum", "detected_cells", "detection_fraction", "pseudobulk_cpm"]
    _columns(genes, [*gene_keys, *gene_metrics], "genes")
    for name in gene_keys:
        genes[name] = _labels(genes[name], f"genes.{name}")
    _require(not genes.duplicated(gene_keys).any(), "genes: duplicate group/mask/sample/gene row")
    _require(set(genes.sample_id) == {reference_sample, query_sample}, "genes: unexpected sample libraries")
    _require(set(genes["mask"]) == set(mask_names), "genes: missing or unexpected masks")
    _require(set(map(tuple, genes[["group_scope", "group_label"]].to_numpy())) == set(gene_groups),
             "genes: required pooled/state groups missing or unknown group labels present")
    selected_genes = genes.loc[genes.gene_symbol.isin(raw_genes)]
    expected_gene_keys = {(scope, label, mask, sample, gene) for scope, label in gene_groups
                         for mask in mask_names for sample in (reference_sample, query_sample) for gene in raw_genes}
    _require(set(map(tuple, selected_genes[gene_keys].to_numpy())) == expected_gene_keys,
             "genes: missing required raw-gene summary rows")
    gene_lookup = selected_genes.set_index(gene_keys)
    gene_row_numbers = {tuple(row[key] for key in gene_keys): index + 1
                        for index, row in enumerate(genes.to_dict("records"))}
    raw_records = []
    for scope, label, mask, sample, gene in sorted(expected_gene_keys):
        key = (scope, label, mask, sample, gene)
        info = selections[key[:4]]
        values = raw[info["selected"], raw_genes.index(gene)]
        n, library, available = info["n_cells"], info["total_selected_library_UMI"], info["mask_available"]
        count = int(values.sum()) if available else None
        detected = int((values > 0).sum()) if available else None
        record = dict(zip(gene_keys, key))
        record["source_row_1based"] = gene_row_numbers[key]
        record.update({name: info[name] for name in info if name != "selected"})
        record.update(raw_count_sum=count, detected_cells=detected,
                      detection_fraction=detected / n if available and n else None,
                      pseudobulk_cpm=count / library * 1e6 if library else None)
        supplied = gene_lookup.loc[key]
        for field in [*gene_metrics, "n_input_cells", "n_metadata_known", "n_metadata_missing", "mask_status"]:
            if field in supplied:
                _equal(supplied[field], record[field], f"genes {key} {field}")
        if "mask_available" in supplied:
            actual = _booleans(pd.Series([supplied.mask_available]), f"genes {key} mask_available")[0]
            _require(actual == available, f"genes {key}: mask_available disagrees with cell membership")
        if "mapping_status" in supplied:
            _equal(supplied.mapping_status, "mapped", f"genes {key} mapping_status")
        raw_records.append(record)

    module_lookup = modules.set_index(module_keys)
    frozen_memberships = {}
    for key, row in module_lookup.iterrows():
        scope, label, mask, module = key
        infos, finite = [], []
        for side, sample in (("reference", reference_sample), ("query", query_sample)):
            info = selections[(scope, label, mask, sample)]
            infos.append(info)
            for field in ("n_cells", "mask_status", "n_metadata_missing"):
                _equal(row[f"{field}_{side}"], info[field], f"modules {key} {field}_{side}")
            n_finite = int(row[f"n_finite_scores_{side}"])
            finite.append(n_finite)
            _require(n_finite <= info["n_cells"], f"modules {key}: finite score count exceeds retained cells")
            median = row[f"median_score_{side}"]
            _require(pd.isna(median) if n_finite == 0 else np.isfinite(float(median)),
                     f"modules {key}: median availability disagrees with finite score count")
        expected = ("unavailable" if not all(info["mask_available"] for info in infos) else
                    "missing" if min(info["n_cells"] for info in infos) == 0 or min(finite) == 0 else
                    "low_cells" if min([*(info["n_cells"] for info in infos), *finite]) < min_cells else "comparable")
        _equal(row.comparison_status, expected, f"modules {key} comparison_status")
        delta = float(row.median_score_query) - float(row.median_score_reference) if expected == "comparable" else None
        _equal(row.median_score_delta, delta, f"modules {key} median_score_delta")
        signature = (module, tuple(sorted(cells.loc[groups[(scope, label)] & masks[mask][0], "cell_id"])))
        if signature in frozen_memberships:
            previous = frozen_memberships[signature]
            for field in [f"{name}_{side}" for name in ("n_cells", "n_finite_scores", "median_score")
                          for side in ("reference", "query")]:
                _equal(row[field], None if pd.isna(previous[field]) else previous[field],
                       f"modules {key}: identical cell membership disagrees in {field}")
        else:
            frozen_memberships[signature] = row

    equivalence = {}
    for mask, (passed, _) in masks.items():
        signature = tuple(sorted(cells.loc[passed, "cell_id"]))
        if signature not in equivalence:
            equivalence[signature] = {"canonical_mask": mask, "masks": [], "membership_sha256": _digest(signature),
                                      "n_cells": len(signature),
                                      "n_cells_reference": int((passed & cells.sample_id.eq(reference_sample)).sum()),
                                      "n_cells_query": int((passed & cells.sample_id.eq(query_sample)).sum())}
        equivalence[signature]["masks"].append(mask)
    sensitivity = []
    module_row_numbers = {tuple(row[key] for key in module_keys): index + 1
                          for index, row in enumerate(modules.to_dict("records"))}
    for (scope, label), members in groups.items():
        for module in module_ids:
            distinct, comparable_masks = {}, []
            for mask in mask_names:
                row = module_lookup.loc[(scope, label, mask, module)]
                signature = tuple(sorted(cells.loc[members & masks[mask][0], "cell_id"]))
                if signature in distinct:
                    first = distinct[signature][1]
                    for field in [f"{name}_{side}" for name in ("n_cells", "n_finite_scores", "median_score")
                                  for side in ("reference", "query")]:
                        _equal(row[field], None if pd.isna(first[field]) else first[field],
                               f"modules {(scope, label, module)}: duplicate-membership mask {mask} {field}")
                    if row.comparison_status == "comparable" or first.comparison_status == "comparable":
                        _equal(row.comparison_status, first.comparison_status, f"modules duplicate masks {mask} status")
                        _equal(row.median_score_delta, float(first.median_score_delta), f"modules duplicate masks {mask} delta")
                else:
                    distinct[signature] = (mask, row)
                if row.comparison_status == "comparable":
                    comparable_masks.append(mask)
            comparable = [(name, float(row.median_score_delta)) for name, row in distinct.values()
                          if row.comparison_status == "comparable"]
            values = [delta for _, delta in comparable]
            primary = module_lookup.loc[(scope, label, "all", module)]
            sensitivity.append({"group_scope": scope, "group_label": label, "module_id": module,
                                "source_row_1based": module_row_numbers[(scope, label, "all", module)],
                                "all_delta": None if pd.isna(primary.median_score_delta) else float(primary.median_score_delta),
                                "n_comparable_masks": len(comparable_masks),
                                "n_distinct_comparable_cell_sets": len(comparable),
                                "canonical_masks": [name for name, _ in comparable],
                                "delta_min": min(values) if values else None,
                                "delta_max": max(values) if values else None,
                                "sign_status": _sign_status(values)})
    raw_lookup = {tuple(record[name] for name in gene_keys): record for record in raw_records}
    contrasts = []
    for scope, label in gene_groups:
        for mask in mask_names:
            for gene in raw_genes:
                reference = raw_lookup[(scope, label, mask, reference_sample, gene)]
                query = raw_lookup[(scope, label, mask, query_sample, gene)]
                eligible = min(reference["n_cells"], query["n_cells"]) >= min_cells
                status = ("unavailable" if not reference["mask_available"] or not query["mask_available"] else
                          "missing" if min(reference["n_cells"], query["n_cells"]) == 0 else
                          "low_cells" if not eligible else "comparable")
                record = {"group_scope": scope, "group_label": label, "mask": mask, "gene_symbol": gene,
                          "reference_sample": reference_sample, "query_sample": query_sample,
                          "comparison_status": status, "n_cells_reference": reference["n_cells"], "n_cells_query": query["n_cells"]}
                for field in ("detection_fraction", "pseudobulk_cpm"):
                    a, b = reference[field], query[field]
                    record[f"{field}_delta"] = b - a if status == "comparable" and a is not None and b is not None else None
                contrasts.append(record)
    return {
        "checks": {"n_cells": len(cells), "n_raw_genes": len(raw_genes), "n_masks": len(mask_names),
                   "n_distinct_mask_cell_sets": len(equivalence), "n_gene_rows_recomputed": len(raw_records),
                   "n_gene_rows_checked": len(selected_genes), "n_module_rows_checked": len(modules),
                   "n_module_ids": len(module_ids), "n_gene_rows_outside_selected_panel": len(genes) - len(selected_genes)},
        "mask_equivalence": list(equivalence.values()), "raw_gene_summary": raw_records,
        "raw_gene_contrasts": contrasts, "module_sensitivity": sensitivity,
        "limitations": [
            "Two libraries do not establish independent biological replication; no P values, FDR, confidence intervals or causal conclusions are calculated.",
            "Only the supplied selected-gene per-cell UMI panel is independently aggregated; other source gene summaries are not recomputed.",
            "CPM uses supplied full-feature per-cell total_counts. The full expression matrix and denominator's agreement with that matrix are not reverified.",
            "Module medians and finite-score counts are external summaries: arithmetic, membership counts, availability and duplicate-mask consistency are checked, but module scores are not recomputed.",
            "Mask membership is supplied by the external audit; this review checks flags and identities, not the original QC classifier or thresholds.",
            "Distinct cell sets remain dependent same-source sensitivity observations; pooled and whole-population subtype summaries are overlapping views.",
        ],
    }
