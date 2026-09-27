"""Core-only checks for versioned inference contracts and audited result summaries."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities import CapabilityRegistry
from eacbp.capabilities.advanced_statistics import PyDESeq2PseudobulkCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.orchestrator.checkpoint import fingerprint
from eacbp.orchestrator.resume import ResumeManager
from eacbp.orchestrator.planning import resolve_task_contract
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import BiologicalDesign, StudyManifest
from eacbp.schemas.task import (
    AssumptionAssessment,
    AssumptionStatus,
    ContrastSpec,
    InferenceContract,
    ScientificResult,
    ScientificResultStatus,
    TaskContract,
    TaskStatus,
    TaskResult,
)
from eacbp.auditor.advanced_statistics import _inference_identity_errors, _scientific_result_errors


def _inference_contract():
    contrast = ContrastSpec(
        kind="categorical",
        vector=("condition", "B", "A"),
        factor="condition",
        tested_level="B",
        reference_level="A",
        label="B versus A",
        effect_definition="PyDESeq2 log2 fold change for B versus A",
    )
    return InferenceContract(
        scientific_question="Does treatment change expression in microglia?",
        estimand=contrast.effect_definition,
        target_population="microglia",
        target_cell_type="microglia",
        feature_scope="all 2 input genes",
        input_artifact_uris=["adata://small/microglia/v1"],
        independent_unit="donor",
        observation_unit="cell",
        paired=False,
        design_formula="~condition",
        counts_source="layers.counts",
        contrast_spec=contrast,
        alpha=0.05,
        confidence_level=0.95,
        fdr_family="deg:microglia:B versus A",
    )


def _scientific_result(contract):
    return ScientificResult(
        status=ScientificResultStatus.ESTIMATED_SUPPORTED,
        summary=(
            "1 of 2 tested features met the prespecified Benjamini-Hochberg "
            f"q < 0.05 rule for {contract.estimand}."
        ),
        n_features_tested=2,
        n_features_supported=1,
        assumptions=[
            AssumptionAssessment(
                name="donor_independence",
                status=AssumptionStatus.UNKNOWN,
                evidence={"independent_unit": "donor", "design_rank": 2},
                reason="Donor independence is not established by metadata or design rank.",
            ),
            AssumptionAssessment(
                name="raw_integer_counts",
                status=AssumptionStatus.PASSED,
                evidence={"counts_source": contract.counts_source},
                reason="Integer raw counts passed input validation.",
            ),
            AssumptionAssessment(
                name="design_full_rank",
                status=AssumptionStatus.PASSED,
                evidence={"design_rank": 2},
                reason="The independently reconstructed design is full rank.",
            ),
        ],
    )


def _result_fixture():
    inference = _inference_contract()
    result = TaskResult(
        task_id="deg_microglia",
        capability="deg",
        method_used="pydeseq2_pseudobulk_v1",
        status=TaskStatus.SUCCESS,
        metrics={
            "inference_contract_id": inference.contract_id,
            "scientific_result_status": "estimated_supported",
            "n_features_tested": 2,
            "n_features_supported": 1,
            "contrast_spec": inference.contrast_spec.as_dict(),
        },
        inference_contract=inference,
        scientific_result=_scientific_result(inference),
    )
    task = TaskContract(
        task_id="deg_microglia",
        capability="deg",
        method="pydeseq2_pseudobulk_v1",
        input_artifacts=list(inference.input_artifact_uris),
        parameters={
            "condition_a": "A",
            "condition_b": "B",
            "donor_col": "donor",
            "paired": False,
            "target_cell_type": "microglia",
            "scientific_question": inference.scientific_question,
            "fdr_family": inference.fdr_family,
        },
        inference_contract=inference,
    )
    table = pd.DataFrame(
        {
            "gene": ["g1", "g2"],
            "p_value": [0.001, 0.7],
            "fdr_q_value": [0.01, 0.7],
            "significant_at_alpha": [True, False],
        }
    )
    return task, result, table, inference


def test_contract_id_is_stable_and_content_sensitive():
    inference = _inference_contract()
    restored = InferenceContract.model_validate(inference.model_dump(mode="json"))
    assert restored.contract_id == inference.contract_id
    changed = inference.model_dump(mode="json")
    changed["fdr_family"] = "another-family"
    with pytest.raises(ValueError, match="contract_id"):
        InferenceContract.model_validate(changed)


def test_resume_signature_includes_new_contract_and_preserves_legacy_shape():
    task, _, _, inference = _result_fixture()
    task.inference_contract = None
    manifest = StudyManifest(
        study_id="small",
        biological_design=BiologicalDesign(species="Mus musculus", tissue="brain"),
    )
    input_hashes = {"adata://small/microglia/v1": "input-hash"}
    environment = {"python": "3.x"}
    signature = ResumeManager.signature(task, input_hashes, manifest, environment)
    task_payload = task.model_dump()
    task_payload.pop("inference_contract", None)
    legacy = fingerprint(
        {
            "task": task_payload,
            "inputs": input_hashes,
            "manifest": manifest.model_dump(),
            "environment": environment,
        }
    )
    assert signature == legacy
    task.inference_contract = inference
    assert ResumeManager.signature(task, input_hashes, manifest, environment) != signature


def test_auditor_recounts_summary_and_checks_each_required_assessment():
    task, result, table, _ = _result_fixture()
    details = SimpleNamespace(design_rank=2)
    params = {"alpha": 0.05}
    assert not _scientific_result_errors(result, [table], params, result.method_used, details)

    tampered = result.model_copy(
        update={"scientific_result": result.scientific_result.model_copy(update={"summary": "Effects are equivalent."})}
    )
    assert any("summary" in error for error in _scientific_result_errors(tampered, [table], params, result.method_used, details))

    tampered = result.model_copy(
        update={"scientific_result": result.scientific_result.model_copy(update={"status": ScientificResultStatus.ESTIMATED_INCONCLUSIVE})}
    )
    assert any("status" in error for error in _scientific_result_errors(tampered, [table], params, result.method_used, details))

    assessments = list(result.scientific_result.assumptions)
    assessments[1] = assessments[1].model_copy(update={"status": AssumptionStatus.UNKNOWN})
    tampered = result.model_copy(
        update={"scientific_result": result.scientific_result.model_copy(update={"assumptions": assessments})}
    )
    assert any("raw counts" in error for error in _scientific_result_errors(tampered, [table], params, result.method_used, details))

    assessments = list(result.scientific_result.assumptions)
    assessments[2] = assessments[2].model_copy(update={"evidence": {"design_rank": 1}})
    tampered = result.model_copy(
        update={"scientific_result": result.scientific_result.model_copy(update={"assumptions": assessments})}
    )
    assert any("design rank" in error for error in _scientific_result_errors(tampered, [table], params, result.method_used, details))

    tampered = result.model_copy(update={"inference_contract_id": "inference-v1:tampered"})
    errors = _inference_identity_errors(task, tampered, None, None, params, None, result.method_used)
    assert any("inference_contract_id" in error for error in errors)


def test_auditor_counts_only_finite_test_statistics():
    task, result, table, _ = _result_fixture()
    table.loc[1, "p_value"] = np.inf
    table.loc[1, "fdr_q_value"] = np.inf
    table.loc[1, "significant_at_alpha"] = False
    science = result.scientific_result.model_copy(
        update={
            "n_features_tested": 1,
            "n_features_supported": 1,
            "summary": (
                "1 of 1 tested features met the prespecified Benjamini-Hochberg "
                f"q < 0.05 rule for {result.inference_contract.estimand}."
            ),
        }
    )
    metrics = dict(result.metrics)
    metrics.update(n_features_tested=1, n_features_supported=1)
    tampered = result.model_copy(update={"scientific_result": science, "metrics": metrics})
    errors = _scientific_result_errors(tampered, [table], {"alpha": 0.05}, result.method_used, SimpleNamespace(design_rank=2))
    assert not errors


def test_underpowered_advanced_statistics_returns_structured_scientific_failure(tmp_path):
    data = SCData(
        X=np.ones((4, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {
                "condition": ["A", "A", "B", "B"],
                "donor": ["d1", "d1", "d2", "d2"],
            }
        ),
        var=pd.DataFrame(index=["g1", "g2"]),
        layers={"counts": np.array([[1, 2], [2, 2], [3, 1], [2, 3]], dtype=np.int64)},
    )
    artifacts = ArtifactRegistry(str(tmp_path / "artifacts"))
    input_uri = "adata://underpowered/input/v1"
    artifacts.register(input_uri, data, ArtifactType.ANNDATA, "underpowered", "input", "input")
    capabilities = CapabilityRegistry()
    capabilities.register(PyDESeq2PseudobulkCapability())
    contract = TaskContract(
        task_id="underpowered_deg",
        capability="deg",
        method="pydeseq2_pseudobulk_v1",
        input_artifacts=[input_uri],
        parameters={"condition_a": "A", "condition_b": "B", "donor_col": "donor", "paired": False},
    )

    class _Router:
        @staticmethod
        def resolve_method(capability, manifest, state):
            return "pydeseq2_pseudobulk_v1"

    manifest = StudyManifest(
        study_id="underpowered",
        biological_design=BiologicalDesign(species="Mus musculus", tissue="brain"),
    )
    resolve_task_contract(contract, manifest, {}, capabilities, artifacts, router=_Router())
    assert contract.inference_contract is not None

    result = capabilities.execute_contract(contract, artifacts)
    assert result.status == TaskStatus.SCIENTIFIC_FAILURE
    assert result.inference_contract is not None
    assert result.inference_contract_id == result.inference_contract.contract_id
    assert result.scientific_result.status == ScientificResultStatus.ASSUMPTIONS_FAILED
    assert result.scientific_result.n_features_tested == 0
    assert result.output_artifacts == []
