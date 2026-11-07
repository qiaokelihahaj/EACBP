"""Reconcile descriptive evidence without inventing samples or denominators."""

import math
import hashlib
import json

import pandas as pd
import pytest

from eacbp.application.descriptive_checks import check_descriptive_tables


REFERENCE = "con"
QUERY = "cKO"
POPULATION = "epithelial_unresolved"


def _tables():
    """Small CSV-shaped tables with equal-sized but different mask memberships."""
    cells = pd.DataFrame([
        ("c1", REFERENCE, "0", 100, 2, 3),
        ("c2", REFERENCE, "0", 200, 0, 1),
        ("c3", REFERENCE, "1", 300, 4, 0),
        ("q1", QUERY, "0", 400, 3, 2),
        ("q2", QUERY, "0", 500, 5, 0),
        ("q3", QUERY, "1", 600, 0, 1),
    ], columns=["cell_id", "sample_id", "cluster", "total_counts", "G1_raw_UMI", "G2_raw_UMI"])
    memberships = {
        "all": {"c1", "c2", "c3", "q1", "q2", "q3"},
        "first": {"c1", "c2", "q1", "q2"},
        "second": {"c1", "c3", "q1", "q3"},
        "clone": {"c1", "c2", "q1", "q2"},
    }
    for mask, members in memberships.items():
        cells[f"mask_{mask}_pass"] = cells.cell_id.isin(members)
        cells[f"mask_{mask}_metadata_known"] = True
    genes, modules = [], []
    groups = [
        ("pooled", "ALL_REVIEWED_EPITHELIAL", cells),
        ("subtype", POPULATION, cells),
        ("subcluster_state", "subcluster_0", cells[cells.cluster == "0"]),
        ("subcluster_state", "subcluster_1", cells[cells.cluster == "1"]),
    ]
    for scope, label, group in groups:
        for mask in memberships:
            selected = group[group[f"mask_{mask}_pass"]]
            common = {"group_scope": scope, "group_label": label, "mask": mask}
            counts = {}
            for sample in (REFERENCE, QUERY):
                sample_cells = selected[selected.sample_id == sample]
                n = len(sample_cells)
                counts[sample] = n
                total = int(sample_cells.total_counts.sum())
                for gene in ("G1", "G2"):
                    values = sample_cells[f"{gene}_raw_UMI"]
                    detected = int((values > 0).sum())
                    raw_sum = int(values.sum())
                    genes.append({
                        **common, "sample_id": sample, "gene_symbol": gene,
                        "n_cells": n, "total_selected_library_UMI": total,
                        "raw_count_sum": raw_sum, "detected_cells": detected,
                        "detection_fraction": detected / n if n else math.nan,
                        "pseudobulk_cpm": raw_sum / total * 1e6 if total else math.nan,
                    })
            n_ref, n_query = counts[REFERENCE], counts[QUERY]
            status = "missing" if min(n_ref, n_query) == 0 else (
                "low_cells" if min(n_ref, n_query) < 2 else "comparable"
            )
            module = {
                **common, "module_id": "response", "reference_sample": REFERENCE,
                "query_sample": QUERY, "comparison_status": status,
                "median_score_delta": 0.5 if status == "comparable" else math.nan,
            }
            for side, n, median in (("reference", n_ref, 0.25), ("query", n_query, 0.75)):
                module.update({
                    f"n_cells_{side}": n, f"n_finite_scores_{side}": n,
                    f"median_score_{side}": median if n else math.nan,
                    f"mask_status_{side}": "available",  # Both libraries exist in every input state, even if filtered to zero.
                    f"n_metadata_missing_{side}": 0,
                })
            modules.append(module)
    return cells, pd.DataFrame(modules), pd.DataFrame(genes)


def _check(tables):
    return check_descriptive_tables(
        *tables, reference_sample=REFERENCE, query_sample=QUERY,
        population=POPULATION, min_cells=2,
    )


def test_reconciles_real_table_shapes_and_full_library_denominator():
    cells, modules, genes = _tables()
    result = _check((cells, modules, genes))
    assert result["checks"]["n_cells"] == 6
    assert result["checks"]["n_raw_genes"] == 2
    assert result["checks"]["n_gene_rows_checked"] == len(genes)
    assert result["checks"]["n_module_rows_checked"] == len(modules)
    summaries = pd.DataFrame(result["raw_gene_summary"])
    row = summaries[
        (summaries.group_scope == "pooled") & (summaries["mask"] == "all")
        & (summaries.sample_id == REFERENCE) & (summaries.gene_symbol == "G1")
    ].iloc[0]
    assert row.raw_count_sum == 6
    assert row.detected_cells == 2
    assert row.detection_fraction == pytest.approx(2 / 3)
    assert row.total_selected_library_UMI == 600
    assert row.pseudobulk_cpm == pytest.approx(10_000)


def test_mask_equivalence_uses_cell_ids_not_equal_counts():
    result = _check(_tables())
    assert result["checks"]["n_masks"] == 4
    assert result["checks"]["n_distinct_mask_cell_sets"] == 3
    classes = {frozenset(row["masks"]) for row in result["mask_equivalence"]}
    assert classes == {frozenset({"all"}), frozenset({"first", "clone"}), frozenset({"second"})}
    rows = result["mask_equivalence"]
    assert len({row["membership_sha256"] for row in rows}) == 3
    assert sorted(row["n_cells"] for row in rows) == [4, 4, 6]


def test_duplicate_masks_do_not_add_independent_module_support():
    result = _check(_tables())
    rows = pd.DataFrame(result["module_sensitivity"])
    row = rows[(rows.group_scope == "pooled") & (rows.module_id == "response")].iloc[0]
    assert row.n_comparable_masks == 4
    assert row.n_distinct_comparable_cell_sets == 3
    assert row.delta_min == row.delta_max == 0.5


@pytest.mark.parametrize("value", [-1, 0.5, float("nan"), float("inf")])
def test_rejects_invalid_raw_umi(value):
    cells, modules, genes = _tables()
    cells["G1_raw_UMI"] = cells.G1_raw_UMI.astype(float)
    cells.loc[0, "G1_raw_UMI"] = value
    with pytest.raises(ValueError):
        _check((cells, modules, genes))


@pytest.mark.parametrize("column,value", [("total_counts", -1), ("total_counts", 3.5), ("total_counts", 1)])
def test_rejects_invalid_or_inconsistent_full_library_umi(column, value):
    cells, modules, genes = _tables()
    cells[column] = cells[column].astype(float)
    cells.loc[0, column] = value
    with pytest.raises(ValueError):
        _check((cells, modules, genes))


def test_rejects_duplicate_cell_identity():
    cells, modules, genes = _tables()
    cells.loc[3, "cell_id"] = cells.loc[0, "cell_id"]
    with pytest.raises(ValueError):
        _check((cells, modules, genes))


@pytest.mark.parametrize("sample", ["", None, "unexpected"])
def test_rejects_missing_or_unrecognized_samples(sample):
    cells, modules, genes = _tables()
    cells.loc[0, "sample_id"] = sample
    with pytest.raises(ValueError):
        _check((cells, modules, genes))


def test_rejects_mask_without_metadata_known_flag():
    cells, modules, genes = _tables()
    cells = cells.drop(columns=["mask_first_metadata_known"])
    with pytest.raises(ValueError):
        _check((cells, modules, genes))


@pytest.mark.parametrize("column,value", [
    ("mask_first_metadata_known", False),
    ("mask_first_metadata_known", None),
    ("mask_first_pass", "perhaps"),
    ("mask_all_pass", False),
])
def test_unknown_or_invalid_selection_cannot_silently_pass(column, value):
    cells, modules, genes = _tables()
    cells[column] = cells[column].astype(object)
    cells.loc[0, column] = value
    with pytest.raises(ValueError):
        _check((cells, modules, genes))


@pytest.mark.parametrize("column", [
    "raw_count_sum", "detected_cells", "detection_fraction",
    "total_selected_library_UMI", "pseudobulk_cpm", "n_cells",
])
def test_rejects_changed_gene_evidence(column):
    cells, modules, genes = _tables()
    genes.loc[0, column] += 1
    with pytest.raises(ValueError):
        _check((cells, modules, genes))


@pytest.mark.parametrize("table_index", [1, 2], ids=["modules", "genes"])
@pytest.mark.parametrize("mutation", ["duplicate", "omit", "omit_all_state_rows"])
def test_required_source_rows_are_unique_and_complete(table_index, mutation):
    tables = list(_tables())
    table = tables[table_index]
    if mutation == "duplicate":
        table = pd.concat([table, table.iloc[[0]]], ignore_index=True)
    elif mutation == "omit":
        table = table.iloc[1:].copy()
    else:
        table = table[table.group_scope != "subcluster_state"].copy()
    tables[table_index] = table
    with pytest.raises(ValueError):
        _check(tables)


@pytest.mark.parametrize("mutation", ["delta", "minimum", "finite_count"])
def test_module_delta_and_minimum_cell_claims_must_agree(mutation):
    cells, modules, genes = _tables()
    if mutation == "delta":
        modules.loc[0, "median_score_delta"] = -0.5
    elif mutation == "minimum":
        index = modules.index[modules.comparison_status == "low_cells"][0]
        modules.loc[index, "comparison_status"] = "comparable"
        modules.loc[index, "median_score_delta"] = 0.5
    else:
        modules.loc[0, "n_finite_scores_reference"] = 4
    with pytest.raises(ValueError):
        _check((cells, modules, genes))


def test_csv_boolean_strings_supported_without_mutating_inputs():
    cells, modules, genes = _tables()
    for column in cells.columns:
        if column.startswith("mask_"):
            cells[column] = cells[column].map({True: "True", False: "False"})
    before = [table.copy(deep=True) for table in (cells, modules, genes)]
    _check((cells, modules, genes))
    for original, table in zip(before, (cells, modules, genes)):
        pd.testing.assert_frame_equal(original, table)


def test_whole_population_alias_cannot_disagree_with_pooled_scores():
    cells, modules, genes = _tables()
    index = modules.index[(modules.group_scope == "subtype") & (modules["mask"] == "all")][0]
    # Preserve arithmetic and cell counts while corrupting the shared score basis.
    modules.loc[index, ["median_score_reference", "median_score_query"]] += 1.0
    with pytest.raises(ValueError, match="identical cell membership"):
        _check((cells, modules, genes))


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(path, value):
    path.write_text(json.dumps(value), encoding="utf-8")
    return _sha(path)


def _bundle(tmp_path):
    cells, modules, genes = _tables()
    files = {}
    for role, frame in (("cells", cells), ("module_summaries", modules), ("gene_summaries", genes)):
        path = tmp_path / (role + ".csv")
        frame.to_csv(path, index=False)
        files[role] = {"path": path.name, "sha256": _sha(path)}
    definitions = {"module_set_id": "test", "version": "1", "modules": [{"id": "response"}]}
    definition_sha = _save(tmp_path / "modules.json", definitions)
    scoring = {"status": "complete", "inferential_statistics": False, "module_file_sha256": definition_sha,
               "module_set": "test", "module_version": "1", "normalization": "joint log1p",
               "seed": 0, "score_parameters": {"ctrl_size": 50}, "score_background_gene_count": 100,
               "n_epithelial_cells": len(cells)}
    scoring_sha = _save(tmp_path / "scoring.json", scoring)
    sensitivity = {"status": "complete", "inferential_statistics": False,
                   "source_scoring_module_sha256": definition_sha, "modules": {"sha256": definition_sha},
                   "source_scoring_manifest": {"sha256": scoring_sha},
                   "parameters": {"reference_sample": REFERENCE, "query_sample": QUERY, "min_cells": 2, "state_column": "cluster"},
                   "original_scores_unchanged": True,
                   "source_score_basis": {key: scoring[key] for key in ("normalization", "seed", "score_parameters", "score_background_gene_count")},
                   "n_source_cells": len(cells), "n_module_delta_rows": len(modules),
                   "score_columns": ["module_response"]}
    _save(tmp_path / "sensitivity.json", sensitivity)
    for role in ("modules", "scoring", "sensitivity"):
        path = tmp_path / (role + ".json")
        files[role] = {"path": path.name, "sha256": _sha(path)}
    spec = {"schema_version": 1, "format": "frozen_scoring_tables_v1", "study_id": "small_export",
            "population": POPULATION, "reference_sample": REFERENCE, "query_sample": QUERY,
            "min_cells": 2, "design": {"observation_unit": "library", "biological_replication": "unverified"},
            "sources": files}
    path = tmp_path / "bundle.json"
    _save(path, spec)
    return path


def test_pinned_bundle_review_and_report_preserve_scientific_boundary(tmp_path):
    from eacbp.application.descriptive_review import review_descriptive_bundle, format_descriptive_review
    bundle = _bundle(tmp_path)
    result = review_descriptive_bundle(bundle, _sha(bundle))
    assert result["review_status"] == "passed"
    assert result["native_evidence_admitted"] is False
    assert result["inference_status"] == "not_estimable"
    assert result["raw_recomputation"] == "selected_exported_genes_only"
    assert len(result["sources"]) == 6
    assert all(row["source_row_1based"] >= 1 for row in result["checks"]["raw_gene_summary"])
    report = format_descriptive_review(result)
    for phrase in ("观察单位为文库", "独立生物学重复", "不输出组间 P 值", "未重新计算模块分数", "状态内比较", "3 / 3"):
        assert phrase in report
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("role", ["cells", "module_summaries", "gene_summaries", "scoring", "sensitivity", "modules", "bundle"])
def test_pinned_bundle_rejects_changed_source_bytes(tmp_path, role):
    from eacbp.application.descriptive_review import review_descriptive_bundle
    bundle = _bundle(tmp_path)
    pin = _sha(bundle)
    target = bundle if role == "bundle" else tmp_path / json.loads(bundle.read_text())["sources"][role]["path"]
    target.write_bytes(target.read_bytes() + b" ")
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        review_descriptive_bundle(bundle, pin)


@pytest.mark.parametrize("mutation", ["scoring_binding", "module_binding", "contrast", "minimum", "score_basis"])
def test_matching_file_hashes_do_not_override_inconsistent_run_bindings(tmp_path, mutation):
    from eacbp.application.descriptive_review import review_descriptive_bundle
    bundle = _bundle(tmp_path)
    spec = json.loads(bundle.read_text())
    path = tmp_path / "sensitivity.json"
    manifest = json.loads(path.read_text())
    if mutation == "scoring_binding":
        manifest["source_scoring_manifest"]["sha256"] = "a" * 64
    elif mutation == "module_binding":
        manifest["source_scoring_module_sha256"] = "a" * 64
    elif mutation == "contrast":
        manifest["parameters"]["reference_sample"] = QUERY
    elif mutation == "minimum":
        manifest["parameters"]["min_cells"] = 1
    else:
        manifest["source_score_basis"]["seed"] = 9
    spec["sources"]["sensitivity"]["sha256"] = _save(path, manifest)
    pin = _save(bundle, spec)
    with pytest.raises(ValueError):
        review_descriptive_bundle(bundle, pin)


def test_output_preflight_preserves_sources_and_existing_results(tmp_path):
    from eacbp.application.descriptive_review import review_descriptive_bundle, publish_descriptive_review
    bundle = _bundle(tmp_path)
    result = review_descriptive_bundle(bundle, _sha(bundle))
    output, report = tmp_path / "review.json", tmp_path / "report.md"
    report.write_text("previous report", encoding="utf-8")
    before = _sha(bundle)
    for bad in (bundle, report):
        with pytest.raises(ValueError, match="new files"):
            publish_descriptive_review(result, output=output, report=bad)
        assert not output.exists()
    assert report.read_text() == "previous report"
    assert _sha(bundle) == before
    with pytest.raises(ValueError, match="distinct"):
        publish_descriptive_review(result, output=output, report=output)
    assert not output.exists()


def test_cli_review_and_failed_pin_does_not_publish(tmp_path, capsys):
    from eacbp.cli import main
    bundle = _bundle(tmp_path)
    output, report = tmp_path / "review.json", tmp_path / "report.md"
    args = ["review-descriptive", "--bundle", str(bundle), "--bundle-sha256", _sha(bundle),
            "--output", str(output), "--report", str(report)]
    assert main(args) == 0
    assert json.loads(output.read_text(encoding="utf-8"))["review_status"] == "passed"
    assert "观察单位为文库" in report.read_text(encoding="utf-8")
    capsys.readouterr()
    failed = tmp_path / "failed.json"
    assert main(["review-descriptive", "--bundle", str(bundle), "--bundle-sha256", "0" * 64,
                 "--output", str(failed)]) == 1
    assert not failed.exists()
    assert "SHA256 mismatch" in capsys.readouterr().err
