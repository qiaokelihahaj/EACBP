"""Guard tests for the read-only PBMC3k public-data acceptance runner."""

from __future__ import annotations

import hashlib
import json
from pathlib import PurePosixPath, PureWindowsPath
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scripts import verify_public_pbmc as public_pbmc
from scripts.fetch_public_validation_data import DATASETS
from scripts.verify_public_pbmc import (
    PBMC_CELL_TYPE_LABEL_MAP,
    _evaluate_cell_type_annotations,
    _evaluate_reference,
    _prepare_working_input,
    _validate_external_work_dir_path,
    evaluate_saved_annotations,
)


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_working_copy_keeps_raw_cells_features_and_counts_without_mutating_source(tmp_path):
    counts = sparse.csr_matrix(
        np.asarray([[2, 0, 1], [0, 3, 0], [1, 1, 0], [0, 2, 4]], dtype=np.int32)
    )
    source = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(index=["bc-a", "bc-b", "bc-c", "bc-d"]),
        var=pd.DataFrame({"gene_symbols": ["CD3D", "IL7R", "MS4A1"]}, index=["id-a", "id-b", "id-c"]),
    )
    source_path = tmp_path / "raw.h5ad"
    prepared_path = tmp_path / "work" / "prepared.h5ad"
    source.write_h5ad(source_path)
    before = _sha256(source_path)

    info = _prepare_working_input(source_path, prepared_path)
    prepared = ad.read_h5ad(prepared_path)

    assert _sha256(source_path) == before
    assert info["counts_source_before_staging"] == "X"
    assert info["raw_cells_and_features_preserved_in_working_copy"] is True
    assert prepared.shape == source.shape
    assert prepared.obs_names.tolist() == source.obs_names.tolist()
    assert prepared.var_names.tolist() == source.var_names.tolist()
    assert sparse.issparse(prepared.layers["counts"])
    assert (prepared.layers["counts"] != source.X).nnz == 0
    assert prepared.var["gene_name"].tolist() == ["CD3D", "IL7R", "MS4A1"]
    assert prepared.obs.columns.tolist() == []


def test_working_copy_rejects_scaled_matrix_when_no_raw_counts_layer_exists(tmp_path):
    scaled = ad.AnnData(
        X=np.asarray([[0.25, 0.75, 0.0], [0.0, 0.1, 0.9], [0.4, 0.6, 0.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["bc-a", "bc-b", "bc-c"]),
        var=pd.DataFrame({"gene_name": ["CD3D", "IL7R", "MS4A1"]}, index=["g1", "g2", "g3"]),
    )
    source_path = tmp_path / "scaled.h5ad"
    prepared_path = tmp_path / "prepared.h5ad"
    scaled.write_h5ad(source_path)

    with pytest.raises(ValueError, match="integer-valued raw counts"):
        _prepare_working_input(source_path, prepared_path)
    assert not prepared_path.exists()


def test_working_copy_uses_and_preserves_an_explicit_raw_counts_layer(tmp_path):
    counts = np.asarray([[1, 2, 0], [0, 1, 3], [2, 0, 1]], dtype=np.int64)
    source = ad.AnnData(
        X=counts.astype(np.float32) / 4,
        obs=pd.DataFrame(index=["bc-a", "bc-b", "bc-c"]),
        var=pd.DataFrame({"gene_name": ["CD3D", "IL7R", "MS4A1"]}, index=["g1", "g2", "g3"]),
    )
    source.layers["counts"] = counts.copy()
    source_path = tmp_path / "with_counts.h5ad"
    prepared_path = tmp_path / "prepared.h5ad"
    source.write_h5ad(source_path)

    info = _prepare_working_input(source_path, prepared_path)
    prepared = ad.read_h5ad(prepared_path)

    assert info["counts_source_before_staging"] == "layers.counts"
    np.testing.assert_array_equal(prepared.layers["counts"], counts)
    np.testing.assert_array_equal(prepared.X, source.X)


def test_gene_name_fallback_uses_var_names_for_marker_matching(tmp_path):
    source = ad.AnnData(
        X=np.asarray([[1, 0, 2], [0, 2, 1], [1, 1, 0]], dtype=np.int32),
        obs=pd.DataFrame(index=["bc-a", "bc-b", "bc-c"]),
        var=pd.DataFrame(index=["CD3D", "IL7R", "MS4A1"]),
    )
    source_path = tmp_path / "var_names.h5ad"
    prepared_path = tmp_path / "prepared.h5ad"
    source.write_h5ad(source_path)

    info = _prepare_working_input(source_path, prepared_path)
    prepared = ad.read_h5ad(prepared_path)

    assert info["gene_name_source"] == "var_names"
    assert prepared.var["gene_name"].tolist() == ["CD3D", "IL7R", "MS4A1"]
    assert set(info["matched_fixed_marker_genes"]) == {"CD3D", "IL7R", "MS4A1"}


def test_reference_metrics_use_only_exact_shared_barcodes_after_clustering():
    clustered = ad.AnnData(
        X=np.zeros((4, 3), dtype=np.float32),
        obs=pd.DataFrame({"cluster": ["0", "1", "1", "0"]}, index=["bc-a", "bc-b", "bc-c", "raw-only"]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    reference = ad.AnnData(
        X=np.zeros((3, 2), dtype=np.float32),
        obs=pd.DataFrame({"louvain": ["A", "B", "A"]}, index=["bc-a", "bc-c", "ref-only"]),
        var=pd.DataFrame(index=["r1", "r2"]),
    )

    metrics = _evaluate_reference(clustered, reference)

    assert metrics["status"] == "computed"
    assert metrics["reference_used_for_pipeline_or_training"] is False
    assert metrics["n_shared_labeled_cells"] == 2
    assert metrics["n_clustered_cells_without_reference_match"] == 2
    assert metrics["n_reference_cells_without_cluster_match"] == 1
    assert metrics["adjusted_rand_index"] == 1.0
    assert metrics["normalized_mutual_information"] == 1.0
    assert metrics["crosstab"] == {
        "reference_labels": ["A", "B"],
        "cluster_labels": ["0", "1"],
        "counts": [[1, 0], [0, 1]],
    }
    assert metrics["reference_coverage"]["status"] == "failed"


def test_reference_coverage_requires_nonmissing_cluster_and_reference_labels():
    clustered = ad.AnnData(
        X=np.zeros((3, 1), dtype=np.float32),
        obs=pd.DataFrame({"cluster": ["0", None, "1"]}, index=["bc-a", "bc-b", "bc-c"]),
    )
    reference = ad.AnnData(
        X=np.zeros((3, 1), dtype=np.float32),
        obs=pd.DataFrame({"louvain": ["B cells", "NK cells", None]}, index=clustered.obs_names),
    )

    metrics = _evaluate_reference(clustered, reference)

    assert metrics["status"] == "computed"
    assert metrics["n_shared_labeled_cells"] == 1
    assert metrics["reference_coverage"] == {
        "status": "failed",
        "criterion": "Every reference barcode has a nonmissing reference label and predicted cluster.",
        "n_reference_cells": 3,
        "n_shared_barcodes": 3,
        "n_shared_labeled_cells": 1,
    }


def test_cell_type_annotation_metrics_use_fixed_eight_class_alias_map():
    clustered = ad.AnnData(
        X=np.zeros((3, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {"cell_type": ["CD4_T_cells", "CD14_monocytes", "Platelets"]},
            index=["bc-a", "bc-b", "bc-c"],
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    reference = ad.AnnData(
        X=np.zeros((3, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"louvain": ["CD4+ T cells", "CD14 Monocytes", "Megakaryocyte"]},
            index=["bc-a", "bc-b", "bc-c"],
        ),
        var=pd.DataFrame(index=["r1"]),
    )

    metrics = _evaluate_cell_type_annotations(clustered, reference)

    assert len(PBMC_CELL_TYPE_LABEL_MAP) == 8
    assert metrics["status"] == "computed"
    assert metrics["reference_used_for_pipeline_or_training"] is False
    assert metrics["n_cells_scored"] == 3
    assert metrics["n_reference_labels_unmapped"] == 0
    assert metrics["n_prediction_labels_unmapped"] == 0
    assert metrics["accuracy"] == 1.0
    assert metrics["macro_f1"] == 1.0
    assert metrics["annotation_ari"] == 1.0
    assert metrics["annotation_crosstab"]["counts"] == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


def test_unknown_pbmc_reference_labels_remain_unmapped_and_are_not_scored():
    clustered = ad.AnnData(
        X=np.zeros((2, 1), dtype=np.float32),
        obs=pd.DataFrame({"cell_type": ["CD4_T_cells", "CD4_T_cells"]}, index=["bc-a", "bc-b"]),
        var=pd.DataFrame(index=["g1"]),
    )
    reference = ad.AnnData(
        X=np.zeros((2, 1), dtype=np.float32),
        obs=pd.DataFrame({"louvain": ["CD4 T cells", "Plasma cells"]}, index=["bc-a", "bc-b"]),
        var=pd.DataFrame(index=["r1"]),
    )

    metrics = _evaluate_cell_type_annotations(clustered, reference)

    assert metrics["n_shared_cells"] == 2
    assert metrics["n_cells_scored"] == 1
    assert metrics["n_reference_labels_unmapped"] == 1
    assert metrics["reference_unmapped_label_counts"] == {"Plasma cells": 1}
    assert "<unmapped reference: Plasma cells>" in metrics["annotation_crosstab"]["reference_labels"]


def test_unknown_prediction_counts_as_error_instead_of_forced_class_mapping():
    clustered = ad.AnnData(
        X=np.zeros((1, 1), dtype=np.float32),
        obs=pd.DataFrame({"cell_type": ["Plasma cells"]}, index=["bc-a"]),
        var=pd.DataFrame(index=["g1"]),
    )
    reference = ad.AnnData(
        X=np.zeros((1, 1), dtype=np.float32),
        obs=pd.DataFrame({"louvain": ["CD4 T cells"]}, index=["bc-a"]),
        var=pd.DataFrame(index=["r1"]),
    )

    metrics = _evaluate_cell_type_annotations(clustered, reference)

    assert metrics["n_cells_scored"] == 1
    assert metrics["n_prediction_labels_unmapped"] == 1
    assert metrics["accuracy"] == 0.0
    assert metrics["macro_f1"] == 0.0
    assert metrics["annotation_crosstab"]["prediction_labels"] == ["<unmapped prediction: Plasma cells>"]


def test_missing_reference_and_prediction_labels_are_not_string_nan_classes():
    clustered = ad.AnnData(
        X=np.zeros((4, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"cell_type": ["CD4_T_cells", "nan", None, "<NA>"]},
            index=["bc-a", "bc-b", "bc-c", "bc-d"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )
    reference = ad.AnnData(
        X=np.zeros((4, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"louvain": [None, "nan", "CD4 T cells", "CD4 T cells"]},
            index=["bc-a", "bc-b", "bc-c", "bc-d"],
        ),
        var=pd.DataFrame(index=["r1"]),
    )

    metrics = _evaluate_cell_type_annotations(clustered, reference)

    assert metrics["n_reference_labels_missing"] == 2
    assert metrics["n_reference_labels_unmapped"] == 0
    assert metrics["n_prediction_labels_missing"] == 3
    assert metrics["n_prediction_labels_unmapped"] == 0
    assert metrics["n_cells_scored"] == 2
    assert metrics["accuracy"] == 0.0
    assert "<missing reference>" in metrics["annotation_crosstab"]["reference_labels"]
    assert "<missing prediction>" in metrics["annotation_crosstab"]["prediction_labels"]
    all_table_labels = (
        metrics["annotation_crosstab"]["reference_labels"]
        + metrics["annotation_crosstab"]["prediction_labels"]
    )
    assert "nan" not in all_table_labels
    assert all("nan" not in value.casefold() for value in all_table_labels)


def test_saved_h5ad_annotation_helper_reads_existing_cluster_output(tmp_path):
    clustered_path = tmp_path / "clustered.h5ad"
    reference_path = tmp_path / "reference.h5ad"
    clustered = ad.AnnData(
        X=np.zeros((2, 1), dtype=np.float32),
        obs=pd.DataFrame({"cell_type": ["B_cells", "NK_cells"]}, index=["bc-a", "bc-b"]),
        var=pd.DataFrame(index=["g1"]),
    )
    reference = ad.AnnData(
        X=np.zeros((2, 1), dtype=np.float32),
        obs=pd.DataFrame({"louvain": ["B cells", "NK cells"]}, index=["bc-a", "bc-b"]),
        var=pd.DataFrame(index=["r1"]),
    )
    clustered.write_h5ad(clustered_path)
    reference.write_h5ad(reference_path)

    metrics = evaluate_saved_annotations(clustered_path, reference_path)

    assert metrics["n_shared_cells"] == 2
    assert metrics["accuracy"] == 1.0
    assert metrics["macro_f1"] == 1.0
    assert metrics["annotation_ari"] == 1.0


def test_public_pbmc_work_dir_path_rules_are_cross_platform_and_filesystem_free():
    windows_repo = PureWindowsPath("C:/workspace/EACBP")
    with pytest.raises(ValueError, match="absolute path"):
        _validate_external_work_dir_path(PureWindowsPath("relative-output"), windows_repo)

    with pytest.raises(ValueError, match="outside the repository/workspace"):
        _validate_external_work_dir_path(windows_repo / "outputs" / "pbmc", windows_repo)

    with pytest.raises(ValueError, match="C: drive"):
        _validate_external_work_dir_path(PureWindowsPath("C:/validation/pbmc"), windows_repo)

    windows_validation = PureWindowsPath("F:/EACBP-public-validation")
    assert _validate_external_work_dir_path(windows_validation, windows_repo) == windows_validation

    posix_repo = PurePosixPath("/workspace/EACBP")
    with pytest.raises(ValueError, match="absolute path"):
        _validate_external_work_dir_path(PurePosixPath("relative-output"), posix_repo)
    with pytest.raises(ValueError, match="outside the repository/workspace"):
        _validate_external_work_dir_path(posix_repo / "outputs" / "pbmc", posix_repo)
    linux_validation = PurePosixPath("/validation")
    assert _validate_external_work_dir_path(linux_validation, PurePosixPath("/workspace")) == linux_validation


def test_public_pbmc_input_checksums_come_from_the_download_catalog():
    catalog = {spec.key: spec.sha256 for spec in DATASETS}
    assert public_pbmc.PUBLIC_INPUT_SHA256["pbmc3k_raw"] == catalog["pbmc3k_raw"]
    assert public_pbmc.PUBLIC_INPUT_SHA256["pbmc3k_processed"] == catalog["pbmc3k_processed"]


@pytest.fixture
def public_pbmc_run(tmp_path, monkeypatch):
    """Exercise the runner and CLI with real H5AD files and a controlled workflow."""
    from eacbp.application import study_service
    from eacbp.artifact import registry as registry_module

    # These small synthetic fixtures use pytest's isolated temporary directory.
    # Production's external-volume policy is exercised separately by the pure
    # path-policy test above; it must not make this CLI contract test depend on
    # an F: drive or fail solely because --basetemp is inside the checkout.
    monkeypatch.setattr(public_pbmc, "_validate_external_work_dir", lambda path: path.resolve())

    source = ad.AnnData(
        X=np.asarray([[2, 0, 1], [0, 3, 0], [1, 1, 2]], dtype=np.int32),
        obs=pd.DataFrame(index=["bc-a", "bc-b", "bc-c"]),
        var=pd.DataFrame(index=["CD3D", "IL7R", "MS4A1"]),
    )
    reference = ad.AnnData(
        X=np.zeros((3, 1), dtype=np.float32),
        obs=pd.DataFrame({"louvain": ["B cells", "NK cells", "CD4 T cells"]}, index=source.obs_names),
    )
    raw_path = tmp_path / "raw.h5ad"
    reference_path = tmp_path / "reference.h5ad"
    source.write_h5ad(raw_path)
    reference.write_h5ad(reference_path)
    # Small fixtures replace the fixed catalog pins only inside this test;
    # the public CLI has no argument or fallback to bypass checksum checks.
    monkeypatch.setattr(public_pbmc, "PUBLIC_INPUT_SHA256", {
        "pbmc3k_raw": _sha256(raw_path),
        "pbmc3k_processed": _sha256(reference_path),
    })
    clustered = source.copy()
    clustered.obs["cluster"] = ["0"] * 3
    clustered.obs["cell_type"] = ["B_cells"] * 3
    context = SimpleNamespace(
        raw_path=raw_path, reference_path=reference_path, work_dir=tmp_path / "run",
        clustered=clustered, calls=[], mutate_reference=False,
    )

    def fake_run_study(**kwargs):
        context.calls.append(kwargs)
        if context.mutate_reference:
            changed_reference = ad.read_h5ad(reference_path)
            changed_reference.obs["louvain"] = ["B cells"] * 3
            changed_reference.write_h5ad(reference_path)
        return {"status": "success"}

    class FakeRegistry:
        def __init__(self, _path):
            pass

        def exists(self, uri):
            return uri.endswith("/annotated/v4")

        def get_scdata(self, _uri):
            return None, SimpleNamespace(to_anndata=lambda: context.clustered)

    monkeypatch.setattr(study_service, "run_study", fake_run_study)
    monkeypatch.setattr(registry_module, "ArtifactRegistry", FakeRegistry)
    return context


def _run_public_pbmc_cli(context):
    exit_code = public_pbmc.main([
        "--data", str(context.raw_path), "--reference", str(context.reference_path),
        "--work-dir", str(context.work_dir),
    ])
    receipt_path, = context.work_dir.glob("pbmc3k-validation-*/public_pbmc_validation.json")
    return exit_code, json.loads(receipt_path.read_text(encoding="utf-8"))


def test_public_pbmc_cli_rejects_replaced_reference_before_workflow(public_pbmc_run):
    context = public_pbmc_run
    context.reference_path.write_bytes(b"replacement reference")

    exit_code, result = _run_public_pbmc_cli(context)

    assert exit_code == 1
    assert result["status"] == "failed"
    assert result["statuses"] == {"engineering": "not_run", "benchmark": "failed"}
    assert "reference SHA256 mismatch" in result["error"]["message"]
    assert result["reference"]["matches_pinned_sha256"] is False
    assert context.calls == []
    assert not list(context.work_dir.glob("*/input/*.h5ad"))


def test_public_pbmc_cli_rejects_replaced_raw_before_workflow(public_pbmc_run):
    context = public_pbmc_run
    context.raw_path.write_bytes(b"replacement raw counts")

    exit_code, result = _run_public_pbmc_cli(context)

    assert exit_code == 1
    assert result["statuses"] == {"engineering": "not_run", "benchmark": "failed"}
    assert "source SHA256 mismatch" in result["error"]["message"]
    assert context.calls == []


def test_public_pbmc_cli_requires_complete_reference_coverage(public_pbmc_run):
    context = public_pbmc_run
    context.clustered = context.clustered[:2, :].copy()

    exit_code, result = _run_public_pbmc_cli(context)

    assert exit_code == 1
    assert result["status"] == "failed"
    assert result["statuses"] == {"engineering": "passed", "benchmark": "failed"}
    assert result["reference_evaluation"]["reference_coverage"]["status"] == "failed"
    assert result["reference_evaluation"]["n_shared_labeled_cells"] == 2
    assert result["reference_evaluation"]["annotation_coverage"]["status"] == "failed"


def test_public_pbmc_cli_requires_all_reference_labels_to_be_scored(public_pbmc_run):
    context = public_pbmc_run
    reference = ad.read_h5ad(context.reference_path)
    reference.obs["louvain"] = ["B cells", "NK cells", "unmapped reference"]
    reference.write_h5ad(context.reference_path)
    public_pbmc.PUBLIC_INPUT_SHA256["pbmc3k_processed"] = _sha256(context.reference_path)

    exit_code, result = _run_public_pbmc_cli(context)

    assert exit_code == 1
    assert result["statuses"] == {"engineering": "passed", "benchmark": "failed"}
    assert result["reference_evaluation"]["reference_coverage"]["status"] == "passed"
    assert result["reference_evaluation"]["annotation_evaluation"]["n_reference_labels_unmapped"] == 1
    assert result["reference_evaluation"]["annotation_coverage"]["status"] == "failed"


def test_public_pbmc_cli_rejects_reference_mutation_during_workflow(public_pbmc_run):
    context = public_pbmc_run
    context.mutate_reference = True

    exit_code, result = _run_public_pbmc_cli(context)

    assert exit_code == 1
    assert result["statuses"] == {"engineering": "passed", "benchmark": "failed"}
    assert result["reference"]["sha256_before"] != result["reference"]["sha256_after"]
    assert result["reference"]["unchanged"] is False
    assert result["integrity_failures"] == ["reference"]


def test_public_pbmc_cli_reports_complete_comparison_without_a_quality_threshold(public_pbmc_run):
    exit_code, result = _run_public_pbmc_cli(public_pbmc_run)

    assert exit_code == 0
    assert result["status"] == "success"
    assert result["statuses"] == {"engineering": "passed", "benchmark": "passed"}
    assert result["source"]["unchanged"] is True
    assert result["reference"]["unchanged"] is True
    assert result["reference_evaluation"]["adjusted_rand_index"] == 0.0
    assert result["reference_evaluation"]["annotation_evaluation"]["accuracy"] < 1.0
    assert result["scientific_interpretation"] == "pending_domain_review"
    assert result["biological_validity_claimed"] is False
