from __future__ import annotations

import math
from pathlib import Path, PurePosixPath, PureWindowsPath

import pytest

from scripts.compare_research_value_reference import (
    _compare_gene_tables,
    _check_primary_fold,
    _prepare_work_dir,
    _validate_external_path,
)


def _result(gene: str, *, effect: object = "1.25") -> dict[str, object]:
    return {
        "gene": gene,
        "log2_fold_change": effect,
        "ci_low": "0.75",
        "ci_high": "1.75",
        "p_value": "0.01",
        "fdr_q_value": "0.02",
    }


def test_external_output_path_validation_is_platform_pure() -> None:
    assert _validate_external_path(
        PurePosixPath("/validation/runs/new"), PurePosixPath("/workspace")
    ) == PurePosixPath("/validation/runs/new")
    assert _validate_external_path(
        PureWindowsPath(r"F:\validation\runs\new"), PureWindowsPath(r"C:\repo")
    ) == PureWindowsPath(r"F:\validation\runs\new")

    with pytest.raises(ValueError, match="absolute"):
        _validate_external_path(PurePosixPath("relative"), PurePosixPath("/workspace"))
    with pytest.raises(ValueError, match="outside"):
        _validate_external_path(PurePosixPath("/workspace/out"), PurePosixPath("/workspace"))
    with pytest.raises(ValueError, match="C:"):
        _validate_external_path(
            PureWindowsPath(r"C:\validation"), PureWindowsPath(r"C:\repo")
        )


def test_work_directory_rejects_relative_path_before_resolving(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="absolute"):
        _prepare_work_dir(Path("relative-run"), project_root=tmp_path)


def test_fold_metadata_rejects_wrong_sample_or_contrast_definition() -> None:
    row = {
        "left_out_donor": "donor-1",
        "fit_status": "estimated",
        "removed_sample_count": "2",
        "design": "~condition + donor",
        "n_pseudobulk_samples": "14",
        "n_donors_condition_a": "7",
        "n_donors_condition_b": "7",
        "paired": "true",
        "condition_a": "stim",
        "condition_b": "ctrl",
    }
    selected, status = _check_primary_fold([row], "donor-1")
    assert selected == [row]
    assert status["status"] == "passed"

    altered = {**row, "n_pseudobulk_samples": "16"}
    _, status = _check_primary_fold([altered], "donor-1")
    assert status["status"] == "failed"
    assert status["fold_field_mismatches"]["n_pseudobulk_samples"]["observed"] == ["16"]


def test_table_comparison_keeps_all_genes_and_matches_nan_and_infinity_masks() -> None:
    primary = [
        _result("G1"),
        _result("G2", effect="NA"),
        _result("G3", effect="inf"),
    ]
    reference = [
        _result("G1", effect=1.25),
        _result("G2", effect=math.nan),
        _result("G3", effect=math.inf),
    ]

    result = _compare_gene_tables(primary, reference, rtol=1e-8, atol=1e-10)

    assert result["status"] == "passed"
    assert result["n_genes_compared"] == 3
    assert result["metrics"]["effect"]["n_nan_primary"] == 1
    assert result["metrics"]["effect"]["n_nan_reference"] == 1
    assert result["metrics"]["effect"]["n_positive_infinity_primary"] == 1
    assert result["metrics"]["effect"]["n_positive_infinity_reference"] == 1
    assert result["metrics"]["effect"]["finite_nan_infinity_masks_equal"] is True


def test_table_comparison_fails_on_mask_or_gene_set_mismatch() -> None:
    primary = [_result("G1"), _result("G2", effect="NA")]
    reference = [_result("G1", effect=1.25), _result("G2", effect=0.5)]

    result = _compare_gene_tables(primary, reference, rtol=1e-8, atol=1e-10)

    assert result["status"] == "failed"
    assert result["metrics"]["effect"]["finite_nan_infinity_masks_equal"] is False
    assert result["metrics"]["effect"]["mask_mismatch_gene_examples"] == ["G2"]

    different_gene_set = _compare_gene_tables(
        [_result("G1")], [_result("G2", effect=1.25)], rtol=1e-8, atol=1e-10
    )
    assert different_gene_set["status"] == "failed"
    assert "no rows were silently dropped" in different_gene_set["reason"]
