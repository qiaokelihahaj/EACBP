import pandas as pd
import pytest
from pathlib import PurePosixPath, PureWindowsPath

from scripts.verify_public_research_value import (
    _candidate_cards,
    _check_fold_row_coverage,
    _known_marker_cards,
    _sha256_assertions,
    _summarize_primary,
    _validate_external_path,
)


def _loo_rows(gene, full_lfc, full_q, folds):
    rows = [{
        "gene": gene,
        "left_out_donor": "__full_model__",
        "fit_status": "full_model",
        "log2_fold_change": full_lfc,
        "ci_low": full_lfc - 0.2,
        "ci_high": full_lfc + 0.2,
        "p_value": 0.001,
        "fdr_q_value": full_q,
    }]
    for donor, lfc, q in folds:
        rows.append({
            "gene": gene,
            "left_out_donor": donor,
            "fit_status": "estimated",
            "log2_fold_change": lfc,
            "ci_low": lfc - 0.3 if pd.notna(lfc) else float("nan"),
            "ci_high": lfc + 0.3 if pd.notna(lfc) else float("nan"),
            "p_value": 0.01 if pd.notna(q) else float("nan"),
            "fdr_q_value": q,
        })
    return rows


def test_primary_summary_surfaces_flip_missing_fold_and_q_change():
    donors = [f"d{i}" for i in range(1, 9)]
    rows = _loo_rows("GENE_FLIP", 2.0, 0.01, [
        ("d1", 1.8, 0.01), ("d2", -0.5, 0.07), ("d3", 2.5, 0.02),
        ("d4", 1.9, 0.01), ("d5", 1.2, 0.03), ("d6", 0.5, 0.04),
        ("d7", 2.1, 0.02),
    ])
    rows += _loo_rows("GENE_ZERO", 1.5, 0.02, [
        ("d1", 0.0, 0.02), ("d2", 1.3, 0.02), ("d3", 1.2, 0.03),
        ("d4", 1.7, 0.04), ("d5", 1.8, 0.01), ("d6", 1.1, 0.02),
        ("d7", 1.5, 0.02), ("d8", 1.6, 0.03),
    ])

    result = _summarize_primary(pd.DataFrame(rows), donors).set_index("gene")
    flip = result.loc["GENE_FLIP"]
    assert flip["finite_folds_of_8"] == "7/8"
    assert flip["same_sign_folds_of_8"] == "6/8"
    assert flip["q_below_alpha_folds_of_8"] == "6/8"
    assert flip["direction_flip_donors"] == ["d2"]
    assert flip["missing_effect_donors"] == ["d8"]
    assert flip["q_status_changed_donors"] == ["d2"]
    assert flip["max_abs_delta_lfc"] == pytest.approx(2.5)
    assert flip["most_influential_donor"] == "d2"
    assert flip["direction_status"] == "direction_flip_and_incomplete_coverage"

    zero = result.loc["GENE_ZERO"]
    assert zero["direction_flip_donors"] == []
    assert zero["zero_effect_donors"] == ["d1"]
    assert zero["direction_status"] == "zero_effect_nonstable"


def test_fit_coverage_allows_nonestimable_gene_values_but_rejects_missing_rows():
    complete = pd.DataFrame([
        {"gene": gene, "left_out_donor": "__full_model__", "fit_status": "full_model", "log2_fold_change": value}
        for gene, value in (("g1", 1.0), ("g2", float("nan")))
    ] + [
        {"gene": gene, "left_out_donor": donor, "fit_status": "estimated", "log2_fold_change": value}
        for donor in ("d1", "d2")
        for gene, value in (("g1", 0.9), ("g2", float("nan")))
    ])
    accepted = _check_fold_row_coverage(complete, ["d1", "d2"], expected_genes=2)
    assert accepted["status"] == "passed"
    assert accepted["non_finite_statistics_are_allowed_when_rows_are_present"] is True

    missing_gene = complete.loc[~((complete["gene"] == "g2") & (complete["left_out_donor"] == "d2"))]
    rejected = _check_fold_row_coverage(missing_gene, ["d1", "d2"], expected_genes=2)
    assert rejected["status"] == "failed"
    assert rejected["donors"]["d2"]["missing_gene_count"] == 1

    missing_fold = complete.loc[complete["left_out_donor"] != "d2"]
    rejected = _check_fold_row_coverage(missing_fold, ["d1", "d2"], expected_genes=2)
    assert rejected["status"] == "failed"
    assert rejected["donors"]["d2"]["complete_gene_table"] is False


def test_exploratory_candidate_cards_exclude_known_positive_benchmark_genes():
    summary = pd.DataFrame([
        {"gene": "ISG15", "full_log2_fold_change": 2.0, "full_fdr_q_value": 0.001,
         "direction_status": "direction_stable_8_of_8", "loo_coverage_status": "complete",
         "n_q_status_changes": 0, "max_abs_delta_lfc": 0.1},
        {"gene": "IFI6", "full_log2_fold_change": 1.4, "full_fdr_q_value": 0.01,
         "direction_status": "direction_flip", "loo_coverage_status": "complete",
         "n_q_status_changes": 1, "max_abs_delta_lfc": 1.5},
        {"gene": "EXPLORATORY1", "full_log2_fold_change": -1.8, "full_fdr_q_value": 0.02,
         "direction_status": "direction_flip", "loo_coverage_status": "complete",
         "n_q_status_changes": 2, "max_abs_delta_lfc": 2.0},
    ])

    cards = _candidate_cards(summary)

    assert cards["gene"].tolist() == ["EXPLORATORY1"]
    assert cards.loc[0, "classification"] == "post_selection_exploratory_sensitivity"
    assert not bool(cards.loc[0, "independent_validation"])
    assert not bool(cards.loc[0, "novelty_claimed"])


def test_known_positive_marker_cards_keep_donor_sensitivity_fields_separate():
    full = pd.DataFrame([{
        "gene": "ISG15", "log2_fold_change": 2.4, "ci_low": 2.0, "ci_high": 2.8,
        "fdr_q_value": 0.001,
    }])
    sensitivity = pd.DataFrame([{
        "gene": "ISG15", "finite_folds_of_8": "8/8", "same_sign_folds_of_8": "8/8",
        "q_below_alpha_folds_of_8": "7/8", "direction_status": "direction_stable_8_of_8",
        "max_abs_delta_lfc": 0.2, "most_influential_donor": "donor3",
        "direction_flip_donors": [], "missing_effect_donors": [], "n_q_status_changes": 1,
    }])
    source_summary = {"marker_direction_check": {
        "status": "passed", "reference": "https://example.test/reference",
        "markers": [{"gene": "ISG15", "status": "positive_direction"}],
    }}

    cards = _known_marker_cards(full, sensitivity, source_summary).set_index("gene")

    assert len(cards) == 7
    assert cards.loc["ISG15", "classification"] == "known_positive_benchmark"
    assert cards.loc["ISG15", "finite_folds_of_8"] == "8/8"
    assert cards.loc["ISG15", "same_sign_folds_of_8"] == "8/8"
    assert cards.loc["ISG15", "q_below_alpha_folds_of_8"] == "7/8"
    assert cards.loc["ISG15", "most_influential_donor"] == "donor3"


def test_provenance_hash_assertions_require_exact_source_files():
    summary = {"files": {
        "data": {"sha256": "data-hash"},
        "manifest": {"sha256": "manifest-file-hash"},
        "table": {"sha256": "table-hash"},
    }}
    checked = _sha256_assertions(
        data_sha256="data-hash", manifest_data_sha256="data-hash",
        manifest_file_sha256="manifest-file-hash",
        source_summary=summary, full_de_sha256="table-hash",
    )
    assert all(checked["checks"].values())

    with pytest.raises(ValueError, match="SHA256 mismatch"):
        _sha256_assertions(
            data_sha256="changed", manifest_data_sha256="data-hash",
            manifest_file_sha256="manifest-file-hash",
            source_summary=summary, full_de_sha256="table-hash",
        )


def test_external_work_directory_policy_is_portable_and_rejects_c_drive_and_repo():
    windows_repo = PureWindowsPath("C:/workspace/EACBP")
    with pytest.raises(ValueError, match="C: drive"):
        _validate_external_path(PureWindowsPath("C:/validation/run"), project_root=windows_repo)
    with pytest.raises(ValueError, match="outside the repository"):
        _validate_external_path(windows_repo / "outputs", project_root=windows_repo)
    external_windows = PureWindowsPath("F:/validation/run")
    assert _validate_external_path(external_windows, project_root=windows_repo) == external_windows
    external_linux = PurePosixPath("/validation/run")
    assert _validate_external_path(external_linux, project_root=PurePosixPath("/workspace/EACBP")) == external_linux
