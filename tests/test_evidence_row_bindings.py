"""Row evidence resolves by scientific identity and fails closed at each boundary."""

import json

import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.auditor.base import ValidationReport
from eacbp.evidence.extraction import (
    EvidenceRowBindingError, _row_digest, extract_evidence,
    validate_evidence_row_binding, validate_evidence_row_structure,
)
from eacbp.evidence.graph import EvidenceGraph
from eacbp.evidence.snapshot import (
    SnapshotIntegrityError, SnapshotValidationError, _stable_fingerprint,
    load_study_snapshot, render_snapshot_report, write_study_snapshot,
)
from eacbp.report.provenance import SentenceProvenanceTracker
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import ClaimNode
from eacbp.schemas.study import BiologicalDesign, StudyManifest
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def _fixture(tmp_path, *, kind="deg", duplicate=False, reverse=False, null_row=True, numeric=False):
    feature = "gene" if kind == "deg" else "source"
    effect = "log2_fold_change" if kind == "deg" else "activity_effect_condition_a_vs_b"
    error = "lfc_se" if kind == "deg" else "activity_se"
    low, high = ("ci_low", "ci_high") if kind == "deg" else ("activity_ci_low", "activity_ci_high")
    rows = [{feature: name, effect: value, error: .2, low: value - .4, high: value + .4,
             "fdr_q_value": q, "condition_a": "stim", "condition_b": "ctrl",
             "contrast_label": "stim versus ctrl", "target_cell_type": "B",
             "status": "estimated", "status_reason": "", "covariates": ""}
            for name, value, q in [("feature1", 1.5, .001), ("feature2", -1.2, .002)]]
    if numeric:
        for row in rows:
            row.update(contrast_kind="numeric", contrast_factor=None,
                       contrast_tested_level=None, contrast_reference_level=None,
                       contrast_label="specified numeric contrast")
    if null_row:
        row = dict(rows[0])
        row.update({feature: "unestimated", effect: float("nan"), error: float("nan"),
                    low: float("nan"), high: float("nan"), "fdr_q_value": float("nan"),
                    "status": "not_estimated", "status_reason": "insufficient_information"})
        rows.append(row)
    if duplicate:
        rows.append(dict(rows[0]))
    if reverse:
        rows.reverse()
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    raw, output = "json://rows/raw/v1", "table://rows/results/v1"
    registry.register(raw, {"input": True}, ArtifactType.JSON, "rows", "ingest", "ingest")
    parameters = {"condition_a": "stim", "condition_b": "ctrl", "target_cell_type": "B", "covariates": []}
    task = TaskContract(task_id="estimate", capability=kind, method="fixture", input_artifacts=[raw],
                        expected_outputs=[output], parameters=parameters)
    transaction = TaskArtifactTransaction(registry)
    transaction.register(output, pd.DataFrame(rows), ArtifactType.TABLE, "rows", task.task_id,
                         "estimate", parameters=parameters, parent_uris=[raw])
    result = TaskResult(task_id=task.task_id, capability=kind, method_used="fixture", status=TaskStatus.SUCCESS,
                        input_artifacts=[raw], output_artifacts=[output],
                        metrics={"is_pseudobulk": True, "covariates": [], "alpha": .05})
    registry.commit_task(transaction, "row_fixture", result)
    transaction.close()
    report = ValidationReport(auditor_name="fixture", target_task_id=result.task_id, overall_passed=True)
    registry.record_audit("row_fixture", task, result, report, auditor_name="fixture")
    manifest = StudyManifest(study_id="rows", data={"raw_artifact_uri": raw},
                             biological_design=BiologicalDesign(species="human", tissue="test"))
    return registry, task, result, report, manifest


def _nodes(fixture):
    registry, task, result, report, _ = fixture
    return extract_evidence(task, result, report, registry)


def _graph(nodes):
    graph = EvidenceGraph()
    for node in nodes:
        graph.add_evidence(node)
    graph.add_claim(ClaimNode(claim_id="C1", statement="An observed result with a contradicting result.",
                             support_evidence_ids=[nodes[0].evidence_id],
                             contradiction_evidence_ids=[nodes[1].evidence_id]))
    return graph


def _write(tmp_path, fixture, nodes):
    registry, _, result, report, manifest = fixture
    path = tmp_path / "study.snapshot.json"
    write_study_snapshot(path, manifest=manifest, config={}, summary={"status": "success"},
                         evidence_graph=_graph(nodes), artifact_registry=registry,
                         task_history=[result], audit_reports=[report])
    return path


@pytest.mark.parametrize("kind", ["deg", "functional_activity"])
def test_real_registered_table_rows_bind_statistics_context_and_null_status(tmp_path, kind):
    fixture = _fixture(tmp_path, kind=kind)
    nodes = _nodes(fixture)
    assert len(nodes) == 2
    for node in nodes:
        selector = node.row_selector
        assert selector.artifact_uri == fixture[2].output_artifacts[0]
        assert selector.key["target_cell_type"] == "B"
        assert node.biological_context["target_cell_type"] == "B"
        assert node.metrics["status_reason"] == "" and node.metrics["covariates"] == []
        assert validate_evidence_row_binding(node, fixture[0], fixture[2]) == "verified"
        assert validate_evidence_row_binding(node, fixture[0]) == "verified"
    # Not-estimated rows retain their NaNs in the actual scientific table;
    # they are not converted into inferential row claims or invented estimates.
    source = fixture[0].get(fixture[2].output_artifacts[0])[1]
    assert source.loc[source.status.eq("not_estimated"), "fdr_q_value"].isna().all()
    path = _write(tmp_path, fixture, nodes)
    restored = load_study_snapshot(path, fixture[0])
    assert all(node.row_selector for node in restored.evidence_graph.evidence_nodes.values())
    markdown = render_snapshot_report(path, fixture[0])
    assert "An observed result" in markdown
    assert "Contradicting Evidence" in markdown
    assert nodes[1].evidence_id in markdown
    assert "feature2" in markdown and "verified" in markdown


@pytest.mark.parametrize("kind", ["deg", "functional_activity"])
def test_row_identity_survives_sorting_without_accepting_changed_artifact_hashes(tmp_path, kind):
    first = _fixture(tmp_path / "first", kind=kind)
    second = _fixture(tmp_path / "second", kind=kind, reverse=True)
    nodes = _nodes(first)
    for node in nodes:
        assert validate_evidence_row_binding(node, second[0], second[2]) == "verified"
    assert first[0].get_metadata(first[2].output_artifacts[0]).sha256_hash != second[0].get_metadata(second[2].output_artifacts[0]).sha256_hash
    path = _write(tmp_path, first, nodes)
    with pytest.raises(SnapshotIntegrityError, match="metadata changed"):
        load_study_snapshot(path, second[0])


@pytest.mark.parametrize("kind", ["deg", "functional_activity"])
def test_duplicate_scientific_key_is_rejected_during_extraction(tmp_path, kind):
    fixture = _fixture(tmp_path, kind=kind, duplicate=True)
    with pytest.raises(EvidenceRowBindingError, match="uniquely.*2 rows"):
        _nodes(fixture)


@pytest.mark.parametrize("name,value", [("gene", "invented"), ("log2_fold_change", 9.0), ("lfc_se", .7),
                                       ("ci_low", -.2), ("ci_high", 8.0), ("fdr_q_value", .04),
                                       ("status", "not_estimated"), ("contrast_label", "ctrl versus stim"),
                                       ("target_cell_type", "T")])
def test_metric_and_identity_tampering_is_rejected(tmp_path, name, value):
    fixture = _fixture(tmp_path)
    node = _nodes(fixture)[0]
    node.metrics[name] = value
    with pytest.raises(EvidenceRowBindingError):
        validate_evidence_row_binding(node, fixture[0], fixture[2])


@pytest.mark.parametrize("name,value", [("gene", "invented"), ("contrast", "ctrl versus stim"),
                                       ("direction", "lower"), ("target_cell_type", "T")])
def test_biological_context_tampering_is_rejected(tmp_path, name, value):
    fixture = _fixture(tmp_path)
    node = _nodes(fixture)[0]
    node.biological_context[name] = value
    with pytest.raises(EvidenceRowBindingError, match="context"):
        validate_evidence_row_binding(node, fixture[0], fixture[2])


def test_recomputed_binding_digest_cannot_replace_independent_source_row_check(tmp_path):
    fixture = _fixture(tmp_path)
    nodes = _nodes(fixture)
    nodes[0].metrics["lfc_se"] = .9
    nodes[0].row_selector.row_sha256 = _row_digest(nodes[0].row_selector, nodes[0].metrics, nodes[0].biological_context)
    assert validate_evidence_row_structure(nodes[0]) == "structure_checked"
    with pytest.raises(EvidenceRowBindingError, match="lfc_se"):
        validate_evidence_row_binding(nodes[0], fixture[0], fixture[2])
    with pytest.raises(SnapshotIntegrityError, match="lfc_se"):
        _write(tmp_path, fixture, nodes)


def test_loaded_and_rendered_snapshots_revalidate_forged_but_self_consistent_row_values(tmp_path):
    fixture = _fixture(tmp_path)
    nodes = _nodes(fixture)
    path = _write(tmp_path, fixture, nodes)
    payload = json.loads(path.read_text(encoding="utf-8"))
    raw = next(item for item in payload["evidence"]["nodes"] if item["evidence_id"] == nodes[0].evidence_id)
    raw["metrics"]["lfc_se"] = .9
    nodes[0].metrics["lfc_se"] = .9
    raw["row_selector"]["row_sha256"] = _row_digest(nodes[0].row_selector, nodes[0].metrics, nodes[0].biological_context)
    payload["fingerprints"]["evidence"] = _stable_fingerprint(payload["evidence"])
    path.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
    # Offline loading validates structure only; a source-table proof needs a registry.
    assert load_study_snapshot(path)
    with pytest.raises(SnapshotIntegrityError, match="lfc_se"):
        load_study_snapshot(path, fixture[0])
    with pytest.raises(SnapshotIntegrityError, match="lfc_se"):
        render_snapshot_report(path, fixture[0])


def test_offline_load_rejects_changed_metrics_even_with_recomputed_section_fingerprint(tmp_path):
    fixture = _fixture(tmp_path)
    path = _write(tmp_path, fixture, _nodes(fixture))
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["evidence"]["nodes"][0]["metrics"]["lfc_se"] = .9
    payload["fingerprints"]["evidence"] = _stable_fingerprint(payload["evidence"])
    path.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
    with pytest.raises(SnapshotValidationError, match="row binding fingerprint"):
        load_study_snapshot(path)


def test_source_scientific_fields_cannot_be_removed_from_binding(tmp_path):
    fixture = _fixture(tmp_path)
    node = _nodes(fixture)[0]
    node.metrics.pop("lfc_se")
    node.row_selector.metric_columns.pop("lfc_se")
    node.row_selector.row_sha256 = _row_digest(node.row_selector, node.metrics, node.biological_context)
    with pytest.raises(EvidenceRowBindingError, match="omits source scientific metric lfc_se"):
        validate_evidence_row_binding(node, fixture[0], fixture[2])


def test_support_and_contradiction_cards_have_equal_source_validation(tmp_path):
    fixture = _fixture(tmp_path)
    nodes = _nodes(fixture)
    tracker = SentenceProvenanceTracker(_graph(nodes), fixture[0], [fixture[2]])
    card = tracker.resolve_claim_provenance_card("C1")
    assert card["evidence_count"] == card["contradiction_evidence_count"] == 1
    assert card["evidence_chain"][0]["row_verification"] == "verified"
    assert card["contradiction_evidence_chain"][0]["row_verification"] == "verified"
    assert card["contradiction_evidence_chain"][0]["row_selector"]["key"]["gene"] == "feature2"
    nodes[1].metrics["lfc_se"] = .9
    nodes[1].row_selector.row_sha256 = _row_digest(nodes[1].row_selector, nodes[1].metrics, nodes[1].biological_context)
    with pytest.raises(EvidenceRowBindingError, match="lfc_se"):
        tracker.resolve_claim_provenance_card("C1")


def test_missing_contradiction_is_not_silently_omitted(tmp_path):
    fixture = _fixture(tmp_path)
    graph = _graph(_nodes(fixture))
    graph.claim_nodes["C1"].contradiction_evidence_ids.append("missing")
    with pytest.raises(ValueError, match="missing evidence"):
        SentenceProvenanceTracker(graph, fixture[0], [fixture[2]]).resolve_claim_provenance_card("C1")


def test_legacy_nodes_without_selector_remain_readable_but_do_not_claim_row_verification(tmp_path):
    fixture = _fixture(tmp_path)
    nodes = _nodes(fixture)
    for node in nodes:
        node.row_selector = None
    path = _write(tmp_path, fixture, nodes)
    restored = load_study_snapshot(path, fixture[0])
    card = SentenceProvenanceTracker(restored.evidence_graph, fixture[0], [fixture[2]]).resolve_claim_provenance_card("C1")
    assert card["evidence_chain"][0]["row_selector"] is None
    assert card["contradiction_evidence_chain"][0]["row_verification"] == "unavailable"


def test_numeric_contrast_nullable_descriptors_are_preserved_without_masking_statistics(tmp_path):
    fixture = _fixture(tmp_path, numeric=True)
    nodes = _nodes(fixture)
    assert nodes[0].metrics["contrast_factor"] is None
    assert nodes[0].metrics["lfc_se"] == .2
    path = _write(tmp_path, fixture, nodes)
    assert load_study_snapshot(path, fixture[0])


@pytest.mark.parametrize("kind", ["deg", "functional_activity"])
def test_native_sentence_cannot_disagree_with_its_valid_source_row(tmp_path, kind):
    fixture = _fixture(tmp_path, kind=kind)
    nodes = _nodes(fixture)
    nodes[0].summary = "Invented_gene changes by 9 for the reversed contrast."
    # The numerical/key binding remains internally consistent. The sentence
    # still must be rejected when independently checked against the source.
    assert validate_evidence_row_structure(nodes[0]) == "structure_checked"
    with pytest.raises(EvidenceRowBindingError, match="summary differs"):
        validate_evidence_row_binding(nodes[0], fixture[0], fixture[2])
    with pytest.raises(SnapshotIntegrityError, match="summary differs"):
        _write(tmp_path, fixture, nodes)


# Regression: real advanced native outputs put population in a nested
# inference contract, not a flat table column or flat computation parameter.
def _nested_target_fixture(tmp_path, *, kind="deg", target="B cells"):
    from eacbp.schemas.task import InferenceContract

    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    raw, output = "json://nested/raw/v1", "table://nested/results/v1"
    registry.register(raw, {"input": True}, ArtifactType.JSON, "nested", "ingest", "ingest")
    definition = "PyDESeq2 log2 fold change for stim versus ctrl" if kind == "deg" else "donor-level inferred activity difference"
    inference = InferenceContract(
        scientific_question="Compare stimulus and control within the stated population",
        estimand=definition, target_population=target or "all input observations", target_cell_type=target,
        feature_scope="all input genes", input_artifact_uris=[raw], independent_unit="donor",
        observation_unit="cell", paired=True, design_formula="~condition + donor", counts_source="layers.counts",
        contrast_spec={"kind": "categorical", "factor": "condition", "tested_level": "stim",
                       "reference_level": "ctrl", "label": "stim versus ctrl",
                       "vector": ["condition", "stim", "ctrl"], "effect_definition": definition},
        fdr_family="deg:observed_population")
    feature = "gene" if kind == "deg" else "source"
    effect = "log2_fold_change" if kind == "deg" else "activity_effect_condition_a_vs_b"
    error = "lfc_se" if kind == "deg" else "activity_se"
    rows = [{feature: name, effect: value, error: .2, "fdr_q_value": q,
             "condition_a": "stim", "condition_b": "ctrl", "contrast_label": "stim versus ctrl",
             "effect_definition": definition, "alpha": .05, "status": "estimated",
             "status_reason": "", "covariates": ""}
            for name, value, q in [("feature1", 1.5, .001), ("feature2", -1.2, .002)]]
    parameters = {"condition_a": "stim", "condition_b": "ctrl", "covariates": [],
                  "inference_contract": inference.model_dump(mode="json"),
                  "inference_contract_id": inference.contract_id}
    task_parameters = {"covariates": [], **({"target_cell_type": target} if target is not None else {})}
    task = TaskContract(task_id="estimate", capability=kind, method="fixture", input_artifacts=[raw],
                        expected_outputs=[output], parameters=task_parameters, inference_contract=inference)
    transaction = TaskArtifactTransaction(registry)
    transaction.register(output, pd.DataFrame(rows), ArtifactType.TABLE, "nested", task.task_id, "estimate",
                         parameters=parameters, parent_uris=[raw],
                         summary_metrics={"inference_contract_id": inference.contract_id})
    result = TaskResult(task_id=task.task_id, capability=kind, method_used="fixture", status=TaskStatus.SUCCESS,
                        input_artifacts=[raw], output_artifacts=[output], inference_contract=inference,
                        metrics={"is_pseudobulk": True, "covariates": [], "alpha": .05,
                                 **({"target_cell_type": target} if target is not None else {})})
    registry.commit_task(transaction, "nested_target_fixture", result)
    transaction.close()
    report = ValidationReport(auditor_name="fixture", target_task_id=result.task_id, overall_passed=True)
    registry.record_audit("nested_target_fixture", task, result, report, auditor_name="fixture")
    manifest = StudyManifest(study_id="nested", data={"raw_artifact_uri": raw},
                             biological_design=BiologicalDesign(species="human", tissue="test"))
    return registry, task, result, report, manifest


def _nested_graph(nodes, fixture):
    graph = _graph(nodes)
    graph.claim_nodes["C1"].inference_contract_ids = [fixture[2].inference_contract_id]
    return graph


def _nested_write(tmp_path, fixture, nodes):
    registry, _, result, report, manifest = fixture
    path = tmp_path / "study.snapshot.json"
    write_study_snapshot(path, manifest=manifest, config={}, summary={"status": "success"},
                         evidence_graph=_nested_graph(nodes, fixture), artifact_registry=registry,
                         task_history=[result], audit_reports=[report])
    return path


@pytest.mark.parametrize("kind", ["deg", "functional_activity"])
def test_native_nested_inference_target_binds_new_row_evidence(tmp_path, kind):
    fixture = _nested_target_fixture(tmp_path, kind=kind)
    nodes = _nodes(fixture)
    for node in nodes:
        assert node.biological_context["target_cell_type"] == "B cells"
        assert "target_cell_type" not in node.row_selector.key
        assert node.row_selector.artifact_parameters["inference_contract"]["target_cell_type"] == "B cells"
        assert validate_evidence_row_binding(node, fixture[0], fixture[2]) == "verified"
    assert load_study_snapshot(_nested_write(tmp_path, fixture, nodes), fixture[0])


@pytest.mark.parametrize("attack", ["change", "remove"])
def test_native_nested_source_target_rejects_context_tampering(tmp_path, attack):
    fixture = _nested_target_fixture(tmp_path)
    node = _nodes(fixture)[0]
    if attack == "change":
        node.biological_context["target_cell_type"] = "T cells"
    else:
        node.biological_context.pop("target_cell_type")
    with pytest.raises(EvidenceRowBindingError, match="target_cell_type"):
        validate_evidence_row_binding(node, fixture[0], fixture[2])


def test_old_v3_selector_without_nested_copy_is_still_source_verified(tmp_path):
    fixture = _nested_target_fixture(tmp_path)
    node = _nodes(fixture)[0]
    node.row_selector.artifact_parameters.pop("inference_contract")
    node.row_selector.row_sha256 = _row_digest(node.row_selector, node.metrics, node.biological_context)
    assert validate_evidence_row_binding(node, fixture[0], fixture[2]) == "verified"
    node.biological_context["target_cell_type"] = "T cells"
    # No frozen target copy and no target row column. The genuine source
    # metadata must supply the independent target check, never node metrics.
    assert validate_evidence_row_structure(node) == "structure_checked"
    with pytest.raises(EvidenceRowBindingError, match="target_cell_type"):
        validate_evidence_row_binding(node, fixture[0], fixture[2])


@pytest.mark.parametrize("side", [0, 1])
def test_support_and_contradiction_both_reject_old_v3_target_tampering(tmp_path, side):
    fixture = _nested_target_fixture(tmp_path)
    nodes = _nodes(fixture)
    for node in nodes:
        node.row_selector.artifact_parameters.pop("inference_contract")
        node.row_selector.row_sha256 = _row_digest(node.row_selector, node.metrics, node.biological_context)
    tracker = SentenceProvenanceTracker(_nested_graph(nodes, fixture), fixture[0], [fixture[2]])
    card = tracker.resolve_claim_provenance_card("C1")
    assert card["evidence_chain"][0]["row_verification"] == "verified"
    assert card["contradiction_evidence_chain"][0]["row_verification"] == "verified"
    nodes[side].biological_context["target_cell_type"] = "T cells"
    with pytest.raises(EvidenceRowBindingError, match="target_cell_type"):
        tracker.resolve_claim_provenance_card("C1")


@pytest.mark.parametrize("source", ["metrics", "inference"])
def test_actual_nested_metadata_target_must_match_task_target(tmp_path, source):
    fixture = _nested_target_fixture(tmp_path)
    node = _nodes(fixture)[0]
    task_copy = fixture[2].model_copy(deep=True)
    if source == "metrics":
        task_copy.metrics["target_cell_type"] = "T cells"
    else:
        task_copy.inference_contract.target_cell_type = "T cells"
    with pytest.raises(EvidenceRowBindingError, match="source task target"):
        validate_evidence_row_binding(node, fixture[0], task_copy)


def test_offline_v3_read_uses_frozen_metadata_target_when_selector_omits_copy(tmp_path):
    fixture = _nested_target_fixture(tmp_path)
    nodes = _nodes(fixture)
    for node in nodes:
        node.row_selector.artifact_parameters.pop("inference_contract")
        node.row_selector.row_sha256 = _row_digest(node.row_selector, node.metrics, node.biological_context)
    path = _nested_write(tmp_path, fixture, nodes)
    payload = json.loads(path.read_text(encoding="utf8"))
    payload["evidence"]["nodes"][0]["biological_context"]["target_cell_type"] = "T cells"
    payload["fingerprints"]["evidence"] = _stable_fingerprint(payload["evidence"])
    path.write_text(json.dumps(payload, allow_nan=False), encoding="utf8")
    with pytest.raises(SnapshotValidationError, match="target_cell_type"):
        load_study_snapshot(path)


def test_old_untargeted_source_cannot_adopt_unbound_metric_and_context_target(tmp_path):
    from eacbp.evidence.extraction import validate_evidence_target_scope

    fixture = _nested_target_fixture(tmp_path, target=None)
    node = _nodes(fixture)[0]
    assert "target_cell_type" not in node.biological_context
    assert validate_evidence_row_binding(node, fixture[0], fixture[2]) == "verified"
    node.metrics["target_cell_type"] = "T cells"
    node.biological_context["target_cell_type"] = "T cells"
    parameters = fixture[0].get_metadata(node.row_selector.artifact_uri).parameters
    with pytest.raises(EvidenceRowBindingError, match="claimed.*absent"):
        validate_evidence_target_scope(node, parameters, fixture[2])
    with pytest.raises(EvidenceRowBindingError, match="claimed.*absent"):
        validate_evidence_row_binding(node, fixture[0], fixture[2])


def test_selector_free_legacy_does_not_gain_false_target_or_row_verification(tmp_path):
    fixture = _nested_target_fixture(tmp_path)
    nodes = _nodes(fixture)
    for node in nodes:
        node.row_selector = None
    assert load_study_snapshot(_nested_write(tmp_path, fixture, nodes), fixture[0])
    card = SentenceProvenanceTracker(_nested_graph(nodes, fixture), fixture[0], [fixture[2]]).resolve_claim_provenance_card("C1")
    assert card["evidence_chain"][0]["row_verification"] == "unavailable"
    assert card["contradiction_evidence_chain"][0]["row_verification"] == "unavailable"
