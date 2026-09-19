"""Regression coverage for explicit FASTQ ingestion and conservative intent parsing."""

from pathlib import Path

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.quantification import FASTQQuantificationCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.orchestrator.intent import IntentParser
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskStatus


def _register_fastq(registry, study_id, payload):
    uri = f"fastq://{study_id}/raw_reads/v1"
    registry.register(
        uri_str=uri,
        payload=payload,
        artifact_type=ArtifactType.FASTQ,
        study_id=study_id,
        created_by_task="task_ingest",
        operation="raw_fastq_manifest_registration",
    )
    return uri


def _contract(uri, *, mode=None, implementation="kb_python_v1"):
    parameters = {"target_gene": "Kat8", "n_cells": 24, "n_genes": 30}
    if mode is not None:
        parameters["mode"] = mode
    return TaskContract(
        task_id="task_quant",
        capability="quantification",
        method=implementation,
        input_artifacts=[uri],
        expected_outputs=[uri.replace("fastq://", "adata://").replace("raw_reads", "raw")],
        parameters=parameters,
    )


def test_real_mode_never_synthesizes_when_reads_are_missing(tmp_path):
    registry = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    uri = _register_fastq(
        registry,
        "missing_reads",
        {
            "samples": {
                "sample_1": {
                    "R1": "/does/not/exist_R1.fastq.gz",
                    "R2": "/does/not/exist_R2.fastq.gz",
                    "metadata": {"donor": "d1", "condition": "control", "batch": "b1"},
                }
            }
        },
    )
    result = FASTQQuantificationCapability().execute(_contract(uri), registry)
    assert result.status == TaskStatus.EXECUTION_FAILURE
    assert result.output_artifacts == []
    assert result.metrics["is_simulated"] is False
    assert "missing" in result.error_message.lower()
    assert not registry.exists("adata://missing_reads/raw/v1")


def test_demo_mode_is_explicit_and_provenance_is_persisted(tmp_path):
    registry = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    uri = _register_fastq(
        registry,
        "demo_reads",
        {
            "samples": {
                "sample_1": {
                    "R1": ["s_L001_R1_001.fastq.gz", "s_L002_R1_001.fastq.gz"],
                    "R2": ["s_L001_R2_001.fastq.gz", "s_L002_R2_001.fastq.gz"],
                    "metadata": {"donor": "d1", "condition": "control", "batch": "b1"},
                }
            }
        },
    )
    result = FASTQQuantificationCapability().execute(_contract(uri, mode="demo"), registry)
    assert result.status == TaskStatus.SUCCESS
    assert result.metrics["is_simulated"] is True
    assert result.metrics["data_origin"] == "synthetic_demo"
    assert result.metrics["lane_pairs_requested"] == 2
    meta, payload = registry.get(result.output_artifacts[0])
    assert meta.summary_metrics["is_simulated"] is True
    assert meta.summary_metrics["data_origin"] == "synthetic_demo"
    uns = payload.uns if hasattr(payload, "uns") else payload["uns"]
    assert uns["is_simulated"] == True
    assert uns["data_origin"] == "synthetic_demo"


def test_starsolo_is_an_explicit_failure(tmp_path):
    registry = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    uri = _register_fastq(registry, "starsolo", {"samples": {}})
    result = FASTQQuantificationCapability(implementation_id="starsolo_v1").execute(
        _contract(uri, mode="demo", implementation="starsolo_v1"), registry
    )
    assert result.status == TaskStatus.METHOD_FAILURE
    assert "demo execution is unavailable" in result.error_message.lower()
    assert result.output_artifacts == []


def test_intent_keeps_adrenal_cortex_and_unknown_counts_unknown():
    manifest = IntentParser.parse_prompt_to_manifest("Analyze human adrenal cortex")
    assert manifest.biological_design.species == "homo_sapiens"
    assert manifest.biological_design.tissue == "adrenal_cortex"
    assert manifest.biological_design.disease is None
    assert manifest.experimental_design.total_samples == 0
    assert manifest.experimental_design.donor_replicates_per_condition == {}


def test_intent_requires_ad_as_a_token():
    manifest = IntentParser.parse_prompt_to_manifest("Analyze an adeno-associated vector in mouse cortex")
    assert manifest.biological_design.disease is None
