"""A reversed contrast keeps its meaning through audit, evidence and reporting."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pydeseq2")

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor import ScientificAuditor
from eacbp.capabilities.advanced_statistics import PyDESeq2PseudobulkCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.evidence.claim import ClaimEngine
from eacbp.evidence.extraction import extract_evidence
from eacbp.evidence.graph import EvidenceGraph
from eacbp.report.markdown_report import ScientificReportGenerator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest
from eacbp.schemas.task import TaskContract


def test_reverse_contrast_nondefault_alpha_through_scientific_audit_and_report(tmp_path):
    rng = np.random.default_rng(1907)
    obs = pd.DataFrame({
        "donor": np.repeat(["d1", "d2", "d3", "d4"], 8),
        "condition": (["A"] * 4 + ["B"] * 4) * 4,
    })
    means = np.full((len(obs), 20), 30)
    means[obs.condition.eq("A"), 0] = 120
    means[obs.condition.eq("B"), 0] = 5
    counts = rng.poisson(means).astype(np.int64)
    data = SCData(counts, obs, pd.DataFrame(index=[f"g{i}" for i in range(20)]), layers={"counts": counts})
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "adata://inference/raw/v1"
    registry.register(uri, data, ArtifactType.ANNDATA, "inference", "input", "input")
    contract = TaskContract(
        task_id="reverse_deg", capability="deg", method="pydeseq2_pseudobulk_v1",
        input_artifacts=[uri], expected_outputs=["table://inference/deg/v1"],
        parameters={"condition_a": "A", "condition_b": "B", "donor_col": "donor",
                    "paired": True, "contrast": ["condition", "B", "A"], "alpha": 0.1,
                    "n_cpus": 1, "independent_filter": False},
        validation_requirements=["advanced_statistics_integrity", "multiple_testing_correction", "pseudoreplication_audit"],
    )
    result = PyDESeq2PseudobulkCapability().execute(contract, registry)
    audit = ScientificAuditor().audit_task(contract, result, registry)
    assert audit.overall_passed, [(check.check_name, check.message) for check in audit.checks if not check.passed]
    evidence = extract_evidence(contract, result, audit, registry)
    gene = next(node for node in evidence if node.biological_context.get("gene") == "g0")
    assert "lower expression for B versus A" in gene.summary
    assert gene.metrics["alpha"] == 0.1
    assert ClaimEngine.has_valid_statistics(gene)
    manifest = StudyManifest(study_id="inference", biological_design=BiologicalDesign(species="human", tissue="test"), data=DataSpec(raw_artifact_uri=uri))
    report = ScientificReportGenerator(manifest, EvidenceGraph(), registry, [result], [audit]).generate_markdown()
    assert "PyDESeq2 log2 fold change for B versus A" in report
    assert "FDR < 0.1:" in report

    result.metrics["alpha"] = 0.05
    assert not ScientificAuditor().audit_task(contract, result, registry).overall_passed
