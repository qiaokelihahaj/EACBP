"""
Unit tests for the Independent Scientific Auditor suite across computational, statistical, spatial, and perturbation checks.
"""

import pytest
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.auditor import ComputationalValidator, StatisticalValidator, BiologicalValidator, ValidationSeverity


def test_computational_validator_catches_nans(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    validator = ComputationalValidator()

    # AnnData with NaN values
    bad_matrix = np.ones((10, 5))
    bad_matrix[2, 3] = np.nan
    data = SCData(
        X=bad_matrix,
        obs=pd.DataFrame({"cell_id": [f"c_{i}" for i in range(10)]}),
        var=pd.DataFrame({"gene_name": [f"g_{i}" for i in range(5)]}),
    )

    reg.register(
        uri_str="adata://AD/corrupted/v1",
        payload=data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD",
        created_by_task="task_corrupted",
        operation="test",
    )

    contract = TaskContract(task_id="task_corrupted", capability="qc")
    result = TaskResult(
        task_id="task_corrupted",
        status=TaskStatus.SUCCESS,
        capability="qc",
        method_used="test",
        output_artifacts=["adata://AD/corrupted/v1"],
    )

    report = validator.audit(contract, result, reg)
    assert not report.overall_passed
    nan_check = next(c for c in report.checks if c.check_name == "expression_finite_values")
    assert not nan_check.passed
    assert nan_check.metrics["nan_count"] == 1


def test_statistical_validator_flags_pseudoreplication(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    validator = StatisticalValidator()

    # Create DEG table executed at cell-level (exploratory)
    deg_cell_level = pd.DataFrame({
        "gene": ["Apoe", "Trem2"],
        "log2_fold_change": [2.5, 1.9],
        "p_value": [1e-12, 1e-8],
        "fdr_q_value": [1e-11, 1e-7],
        "pseudoreplication_warning": [True, True],
    })

    reg.register(
        uri_str="table://AD/deg_cell_level/v1",
        payload=deg_cell_level,
        artifact_type=ArtifactType.TABLE,
        study_id="AD",
        created_by_task="task_deg_cell",
        operation="differential_expression",
        summary_metrics={"statistical_unit": "single_cell"},
    )

    contract = TaskContract(task_id="task_deg_cell", capability="deg")
    result = TaskResult(
        task_id="task_deg_cell",
        status=TaskStatus.SUCCESS,
        capability="deg",
        method_used="cell_level_mannwhitney",
        output_artifacts=["table://AD/deg_cell_level/v1"],
    )

    report = validator.audit(contract, result, reg)
    pseudo_check = next(c for c in report.checks if c.check_name == "pseudoreplication_audit")
    assert not pseudo_check.passed
    assert pseudo_check.severity == ValidationSeverity.WARNING
    assert pseudo_check.metrics["confirmatory_allowed"] is False


def test_statistical_validator_morans_i_and_gearys_c_bounds(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    validator = StatisticalValidator()

    # Create invalid spatial DEG table with out-of-bounds Moran's I (> 1.0) and negative Geary's C
    invalid_spatial_deg = pd.DataFrame({
        "gene": ["Trem2", "Apoe"],
        "morans_i": [1.5, -0.4],  # 1.5 is > 1.0 (out of bounds)
        "gearys_c": [-0.2, 0.8],  # -0.2 is < 0 (out of bounds)
        "fdr_q_value": [0.001, 0.02],
    })

    reg.register(
        uri_str="table://AD/spatial_deg_invalid/v1",
        payload=invalid_spatial_deg,
        artifact_type=ArtifactType.TABLE,
        study_id="AD",
        created_by_task="task_sp_deg",
        operation="spatial_deg",
    )

    contract = TaskContract(task_id="task_sp_deg", capability="spatial_deg")
    result = TaskResult(
        task_id="task_sp_deg",
        status=TaskStatus.SUCCESS,
        capability="spatial_deg",
        method_used="spatial_deg_morans_i_v1",
        output_artifacts=["table://AD/spatial_deg_invalid/v1"],
    )

    report = validator.audit(contract, result, reg)
    assert not report.overall_passed
    moran_check = next(c for c in report.checks if c.check_name == "morans_i_bounds_check")
    assert not moran_check.passed

    geary_check = next(c for c in report.checks if c.check_name == "gearys_c_bounds_check")
    assert not geary_check.passed


def test_statistical_validator_perturbation_bounds(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    validator = StatisticalValidator()

    # Valid perturbation AnnData
    valid_data = SCData(
        X=np.abs(np.random.rand(10, 5)),
        obs=pd.DataFrame({"cell_id": [f"c_{i}" for i in range(10)]}),
        var=pd.DataFrame({"gene_name": [f"g_{i}" for i in range(5)]}),
    )

    reg.register(
        uri_str="adata://AD/perturb_valid/v1",
        payload=valid_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD",
        created_by_task="task_perturb",
        operation="simulate_genetic_perturbation",
    )

    contract = TaskContract(
        task_id="task_perturb",
        capability="genetic_perturbation_simulation",
        parameters={"network_attenuation": 0.35},
    )
    result = TaskResult(
        task_id="task_perturb",
        status=TaskStatus.SUCCESS,
        capability="genetic_perturbation_simulation",
        method_used="in_silico_crispr_ko_v1",
        output_artifacts=["adata://AD/perturb_valid/v1"],
        metrics={"network_attenuation": 0.35},
    )

    report = validator.audit(contract, result, reg)
    assert report.overall_passed
    shift_check = next(c for c in report.checks if c.check_name == "perturbation_shift_bounds_check")
    assert shift_check.passed
    atten_check = next(c for c in report.checks if c.check_name == "network_attenuation_bounds_check")
    assert atten_check.passed


def test_statistical_validator_epistemic_tagging(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    validator = StatisticalValidator()

    # Untagged report under prior-guided contract
    reg.register(
        uri_str="json://AD/knowledge_untagged/v1",
        payload={"summary": "Unbiased finding without epistemic disclaimer"},
        artifact_type=ArtifactType.JSON,
        study_id="AD",
        created_by_task="task_know",
        operation="knowledge_report",
    )

    contract = TaskContract(
        task_id="task_know",
        capability="knowledge_retrieval",
        parameters={"prior_guided": True, "hypotheses": ["DAM hypothesis"]},
    )
    result = TaskResult(
        task_id="task_know",
        status=TaskStatus.SUCCESS,
        capability="knowledge_retrieval",
        method_used="knowledge_engine_prior_v1",
        output_artifacts=["json://AD/knowledge_untagged/v1"],
        metrics={"summary": "Unbiased finding without epistemic disclaimer"},
    )

    report = validator.audit(contract, result, reg)
    tag_check = next(c for c in report.checks if c.check_name == "epistemic_tagging_check")
    assert not tag_check.passed
    assert tag_check.severity == ValidationSeverity.ERROR
