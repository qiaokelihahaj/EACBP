"""Real file pins and synthetic independent receipts; no scientific data parser."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from eacbp.application.verification import verify_acceptance_index


BRANCHES = ("primary", "resolution_0_5", "resolution_2_0", "clustering_seed_17", "uniform_sample_80pct")
PRIMARY, SECONDARY = "inferred cell type", "cell type"
CAPABILITIES = ("dataset_audit", "qc_review", "qc", "normalization", "integration", "clustering", "cluster_markers")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def disposition(total, kept, raw=None):
    raw = total if raw is None else raw
    return dict(full_source_cells=total, branch_raw_cells=raw, retained_cells=kept,
                omitted_by_frozen_source_selection=total - raw, actual_qc_candidate_cells=raw - kept,
                full_source_retention=kept / total if total else None,
                branch_raw_coverage=raw / total if total else None,
                qc_retention_among_branch_raw=kept / raw if raw else None)


def science_field(ari, *, secondary=False):
    valid, missing, common, classes = (10, 0, 9, 1) if secondary else (8, 2, 8, 3)
    coverage = dict(raw_source_n=10, valid_label_n=valid, missing_label_n=missing,
                    full_source=disposition(10, 9), valid_reference=disposition(valid, common),
                    missing_reference=disposition(missing, 9 - common),
                    per_class={"single": disposition(10, 9)} if secondary else {
                        "A": disposition(6, 6), "unknown": disposition(1, 1), "unassigned": disposition(1, 1)},
                    all_valid_classes_kept_in_denominator=True)
    estimable = not secondary and ari is not None
    agreement = dict(status="estimable" if estimable else "not_estimable", common_valid_label_cells=common,
                     reference_common_classes=classes, prediction_common_classes=2 if estimable or secondary else 1,
                     adjusted_rand_index=ari if estimable else None, normalized_mutual_information=.7 if estimable else None,
                     nmi_average_method="arithmetic")
    values = (True, True, ari >= .5 if estimable else None, True if estimable else None)
    checks = [dict(name=name, passed=value) for name, value in zip(
        ("overall_valid_reference_retention", "every_valid_reference_class_retention", "adjusted_rand_index", "normalized_mutual_information"), values)]
    return dict(status="not_estimable" if not estimable else "passed" if all(values) else "failed",
                coverage=coverage, agreement=agreement, checks=checks,
                failed_checks=[row["name"] for row in checks if row["passed"] is False],
                not_estimable_checks=[row["name"] for row in checks if row["passed"] is None],
                fixed_thresholds=dict(adjusted_rand_index_minimum=.5, normalized_mutual_information_minimum=.6, retention_minimum=.8))


def hypothesis(first, target):
    second = target["agreement"]["adjusted_rand_index"]
    delta = second - first if first is not None and second is not None else None
    improvement = "not_estimable" if delta is None else "passed" if delta >= .05 else "failed"
    combined = "not_estimable" if delta is None or target["status"] == "not_estimable" else "passed" if improvement == target["status"] == "passed" else "failed"
    return dict(hypothesis="Frozen combined hypothesis", reference_field=PRIMARY, primary_ari=first,
                resolution_0_5_ari=second, ari_improvement=delta, minimum_improvement=.05,
                improvement_check=dict(status=improvement, passed=None if delta is None else delta >= .05,
                                       ari_improvement=delta, minimum_improvement=.05),
                target_acceptance=dict(branch_id="resolution_0_5", reference_field=PRIMARY, status=target["status"],
                                       checks=target["checks"], failed_checks=target["failed_checks"],
                                       not_estimable_checks=target["not_estimable_checks"], thresholds=target["fixed_thresholds"]),
                status=combined)


class Fixture:
    def __init__(self, directory):
        self.directory = directory
        self.pins = {}
        self.source = directory / "workspace/source.py"
        self.source_pin = self.write("workspace/source.py", b"old frozen source\n")
        self.data = self.write("data/counts.bin", b"public synthetic counts")
        self.reference = self.write("reference/labels.tsv", b"REFERENCE VALUES MUST NEVER BE PARSED\n")
        old_code = self.write("old_checker.py", b"old independent checker\n")
        native_code = self.write("native_v2.py", b"actual v2 native checker\n")
        external_code = self.write("external_v2.py", b"actual v2 scientific checker\n")
        self.native_code = native_code
        failure = self.write_json("failed_gate.json", dict(status="failed_native_scope", reference_labels_opened=False, external_scores_opened=False))
        diagnosis = self.write_json("diagnosis.json", dict(status="diagnosed_read_only"))
        tests = self.write_json("selftests.json", dict(status="passed"))
        masks = self.write_json("masks.json", dict(reference_values_used_for_masks=False))
        checker_map = {old_code["path"]: old_code["sha256"]}
        self.protocol = dict(status="prespecified_before_native_and_external_scores", new_external_check=dict(project="segerstolpe", reference=self.reference,
            primary_field=PRIMARY, secondary_field=SECONDARY, missing_tokens_casefold=["", "na", "nan", "not available", "not applicable"],
            hypothesis_branch="resolution_0_5", baseline_branch="primary", minimum_ari_improvement=.05,
            minimum_absolute_ari=.5, minimum_nmi=.6, minimum_overall_labelled_source_retention=.8, minimum_per_class_source_retention=.8),
            checker_files=checker_map)
        protocol_pin = self.write_json("protocol.json", self.protocol)
        datasets = []
        for project in ("pbmc3k", "segerstolpe"):
            bindings = self.write_json(project + "/bindings.json", dict(branches=[dict(branch_id=value) for value in BRANCHES]))
            datasets.append(dict(project=project, prepared=dict(bindings_path=bindings["path"], bindings_sha256=bindings["sha256"]),
                                 expected_branch_count=5, original_source=self.data["path"], original_source_sha256=self.data["sha256"]))
        supplement = dict(status="frozen_before_first_native_execution", reference_labels_supplied=False,
                          protocol_path=protocol_pin["path"], protocol_sha256=protocol_pin["sha256"],
                          masks_path=masks["path"], masks_sha256=masks["sha256"], datasets=datasets,
                          checker_files=checker_map, environment={"python": "synthetic"})
        supplement_pin = self.write_json("supplement.json", supplement)
        original = {p["path"]: p["sha256"] for p in (protocol_pin, supplement_pin, masks, failure, diagnosis, old_code)}
        self.addendum = dict(status="validation_only_addendum_post_native_pre_external_scores",
            science_hypothesis_and_native_configuration_unchanged=True, native_rerun_count=0,
            sdrf_reference_values_read_by_score_evaluator=False, original_files=original,
            new_files={p["path"]: p["sha256"] for p in (native_code, external_code)},
            original_failed_gate=failure, diagnosis=diagnosis, verification_receipts=[tests])
        addendum_pin = self.write_json("addendum.json", self.addendum)
        latest = []
        for project in ("pbmc3k", "segerstolpe"):
            comparison = self.write_json(project + "/current/comparison.json", dict(status="completed"))
            report = self.write(project + "/current/report.md", b"completed current report\n")
            pointer = dict(comparison=comparison["path"], comparison_sha256=comparison["sha256"],
                           report=report["path"], report_sha256=report["sha256"], invocation_id=project + "-current")
            pointer_pin = self.write_json(project + "/latest_comparison.json", pointer)
            latest.append(dict(project=project, pointer_path=pointer_pin["path"], pointer_sha256=pointer_pin["sha256"],
                               comparison_path=pointer["comparison"], report_path=pointer["report"],
                               comparison_sha256=pointer["comparison_sha256"], report_sha256=pointer["report_sha256"],
                               invocation_id=pointer["invocation_id"], historical_root_comparison_used=False))
        scopes = [dict(project=project, branch_id=branch, run_dir=str(directory / project / "runs" / branch), native_status="success", question_status="complete",
                       native_payload_values_loaded_only_through_get_audited=True,
                       durable_receipts=[dict(task_id="task" + str(i), capability=capability, signature=hashlib.sha256(capability.encode()).hexdigest(),
                                              artifact_hashes={"artifact://synthetic/" + str(i): "sha256:" + self.data["sha256"]}, get_audited_verified=True)
                                         for i, capability in enumerate(CAPABILITIES)])
                  for project in ("pbmc3k", "segerstolpe") for branch in BRANCHES]
        self.native = dict(status="passed", all_10_native_scopes_admitted_and_question_complete=True,
            all_8_pairs_independently_recomputed=True, reference_labels_opened=False, external_scores_opened=False,
            supplement_sha256=supplement_pin["sha256"], protocol_sha256=protocol_pin["sha256"], masks_sha256=masks["sha256"],
            validation_addendum_sha256=addendum_pin["sha256"], validation_addendum_path=addendum_pin["path"],
            protected_files={p["path"]: p["sha256"] for p in (protocol_pin, supplement_pin, masks)},
            code_identity=dict(files={native_code["path"]: native_code["sha256"], self.source_pin["path"]: self.source_pin["sha256"]}),
            native_scopes=scopes, pairs=[dict(project=project, branch_id=branch) for project in ("pbmc3k", "segerstolpe") for branch in BRANCHES[1:]],
            latest_comparisons=latest)
        native_pin = self.write_json("native.json", self.native)
        rows = [dict(project="segerstolpe", branch_id=branch, run_dir=str(directory / "segerstolpe" / "runs" / branch),
                     native_receipts=[dict(task_id="task" + str(i), capability=capability, signature=hashlib.sha256(capability.encode()).hexdigest(),
                                          artifact_hashes={"artifact://synthetic/" + str(i): "sha256:" + self.data["sha256"]}, all_siblings_get_audited=True) for i, capability in enumerate(CAPABILITIES)],
                     reference_fields={PRIMARY: science_field(.4 if branch == "resolution_0_5" else .3), SECONDARY: science_field(None, secondary=True)}) for branch in BRANCHES]
        self.scientific = dict(status="evaluated", engineering_status="passed", reference_labels_opened_after_all_native_scopes_and_pairs=True,
            gate=dict(**native_pin, validation_addendum_sha256=addendum_pin["sha256"], validation_addendum_path=addendum_pin["path"],
                protocol_sha256=protocol_pin["sha256"], supplement_sha256=supplement_pin["sha256"], masks_sha256=masks["sha256"],
                actual_native_verifier_path=native_code["path"], environment=supplement["environment"],
                all_10_native_scopes_currently_get_audited=True, all_8_pairs_independently_recomputed=True),
            protected_files=self.native["protected_files"], code_identity=dict(files={external_code["path"]: external_code["sha256"], self.source_pin["path"]: self.source_pin["sha256"]}),
            branches=rows, primary_hypothesis=hypothesis(.3, rows[1]["reference_fields"][PRIMARY]),
            reference=dict(canonical_path=self.reference["path"], canonical_sha256=self.reference["sha256"],
                           fields=dict(primary=PRIMARY, secondary_descriptive=SECONDARY),
                           missing_tokens_stripped_casefold=self.protocol["new_external_check"]["missing_tokens_casefold"]),
            interpretation="Deposited expression-derived labels; independent biological truth was not established.",
            prior_knowledge_and_exposure={"fully_blinded": False},
            winner_selected=False, primary_replaced=False, default_changed=False)
        science_pin = self.write_json("scientific.json", self.scientific)
        regression = self.write_json("regression.json", dict(status="passed", both_pbmc_primary_native_scopes_and_durable_audits_verified=True,
            code_identity=dict(source_files={old_code["path"]: old_code["sha256"]}), protected_old_files={self.data["path"]: self.data["sha256"]},
            protected_new_files={self.data["path"]: self.data["sha256"]},
            comparison_policy=dict(required_capabilities=list(CAPABILITIES[2:]), raw_count_layers_and_dtype_exact=True,
                scientific_arrays_and_tables_exact=True, cell_and_feature_order_exact=True,
                method_operation_effective_scientific_parameters_must_match=True,
                numeric_rtol=0, numeric_atol=0, embedding_rtol=0, embedding_atol=0),
            projects=[dict(project="pbmc3k", status="passed", failed_capabilities=[], capabilities=[
                dict(capability=name, old_method="same", new_method="same", passed=True, outputs=[
                    dict(passed=True, effective_scientific_parameters_equal=True, old_effective_parameters={}, new_effective_parameters={},
                         semantic=dict(passed=True, X=dict(passed=True, exact=True, rtol=0, atol=0)))]) for name in CAPABILITIES[2:]])]))
        recovery = self.write_json("p4-independent-recovery-v1.json", dict(status="passed_independent_recovery",
            all_expected_fault_and_recovery_phases_verified=True, native_execution_attempts_total=10,
            references_opened=False, projects={"pbmc3k": {"runs": 5}, "segerstolpe": {"runs": 5}},
            supplement=supplement_pin, checker=old_code))
        lanes = []
        for name in ("core", "standard", "advanced-statistics"):
            junit = self.write("ci-" + name + ".xml", b'<testsuites><testsuite tests="2" failures="0" errors="0" skipped="0"/></testsuites>')
            log = self.write("ci-" + name + ".log", b"2 passed\n")
            receipt = self.write_json("ci-" + name + ".json", dict(returncode=0, tests=2, failures=0, errors=0, skipped=0,
                explicit_deselect=[], environment_before=supplement["environment"], environment_after=supplement["environment"],
                junit_path=junit["path"], log_path=log["path"]))
            lanes.append(dict(**receipt, lane=name, tests=2, explicit_deselect=[], junit=junit, log=log))
        self.excluded_log = directory / "publication.log"
        self.excluded_log.write_text("still writing\n")
        self.index = dict(schema_version=1, engineering_status="passed", scientific_primary_hypothesis=deepcopy(self.scientific["primary_hypothesis"]),
            native_scopes=deepcopy(self.native["native_scopes"]), independent_pairs=deepcopy(self.native["pairs"]),
            external_science=deepcopy(rows), external_reference=deepcopy(self.scientific["reference"]),
            independent_native_gate=native_pin, external_science_gate=science_pin, post_native_validation_addendum=addendum_pin,
            preserved_failed_native_gate=failure, primary_scientific_regression=regression,
            fault_recovery={"pbmc3k": {"runs": 5}, "segerstolpe": {"runs": 5}}, ci=dict(total_passed=6, lanes=lanes),
            source_files={self.source_pin["path"]: self.source_pin["sha256"]}, prior_knowledge_and_exposure={"fully_blinded": False},
            scientific_interpretation=self.scientific["interpretation"], biological_validation_established=False,
            publication_logs_outside_index=[str(self.excluded_log)], winner_selected=False, primary_replaced=False, default_changed=False)
        self.index_path = directory / "acceptance.json"
        self.publish()

    def write(self, name, data):
        path = self.directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        pin = dict(path=str(path), sha256=digest(path), bytes=len(data))
        self.pins[str(path)] = pin
        return pin

    def write_json(self, name, value):
        return self.write(name, json.dumps(value, sort_keys=True, allow_nan=False).encode())

    def publish(self):
        self.index["evidence_files"] = list(self.pins.values())
        self.index_path.write_text(json.dumps(self.index, sort_keys=True, allow_nan=False))
        self.index_sha = digest(self.index_path)

    def update_science(self):
        self.index["external_science_gate"] = self.write_json("scientific.json", self.scientific)
        self.index["external_science"] = deepcopy(self.scientific["branches"])
        self.index["scientific_primary_hypothesis"] = deepcopy(self.scientific["primary_hypothesis"])
        self.publish()

    def update_native(self):
        pin = self.write_json("native.json", self.native)
        self.index["independent_native_gate"] = pin
        self.index["native_scopes"] = deepcopy(self.native["native_scopes"])
        self.index["independent_pairs"] = deepcopy(self.native["pairs"])
        self.scientific["gate"].update(pin)
        self.update_science()

    def archive(self):
        target = self.directory / "archive/source.py"
        target.parent.mkdir()
        target.write_bytes(self.source.read_bytes())
        receipt = self.directory / "archive.json"
        receipt.write_text(json.dumps(dict(status="passed", p4_index=dict(path=str(self.index_path), sha256=self.index_sha),
            sources=[dict(original_path=str(self.source), archive_path=str(target), sha256=self.source_pin["sha256"])])))
        return receipt, digest(receipt), target

    def verify(self, **options):
        return verify_acceptance_index(self.index_path, expected_index_sha256=self.index_sha, **options)


@pytest.fixture
def evidence(tmp_path):
    return Fixture(tmp_path)


def test_failed_hypothesis_is_not_an_engineering_or_closure_failure(evidence):
    result = evidence.verify()
    assert result.verification_status == result.engineering.status == "passed"
    assert result.science.status == "failed"
    assert result.science.improvement["status"] == "passed"
    assert result.science.absolute_acceptance["status"] == "failed"
    assert result.science.failure_conditions[0]["condition"] == "adjusted_rand_index"
    assert "minimum=0.5" in result.summary
    assert result.engineering.durable_task_count == 70
    assert not result.labels_read and not result.scores_recomputed and not result.scientific_result_reestablished
    assert result.model_dump(mode="json")["science"]["branches"][0]["reference_fields"][SECONDARY]["agreement"]["adjusted_rand_index"] is None


def test_no_reference_parser_or_business_score_function_is_used(evidence, monkeypatch):
    import pandas
    import sklearn.metrics

    def forbidden(*args, **kwargs):
        raise AssertionError("Scientific values must not be parsed or recomputed")

    monkeypatch.setattr(pandas, "read_csv", forbidden)
    monkeypatch.setattr(sklearn.metrics, "adjusted_rand_score", forbidden)
    monkeypatch.setattr(sklearn.metrics, "normalized_mutual_info_score", forbidden)
    assert evidence.verify().verification_status == "passed"


def test_not_estimable_hypothesis_preserves_nulls_without_engineering_failure(evidence):
    for row in evidence.scientific["branches"]:
        row["reference_fields"][PRIMARY] = science_field(None)
    evidence.scientific["primary_hypothesis"] = hypothesis(None, evidence.scientific["branches"][1]["reference_fields"][PRIMARY])
    evidence.update_science()
    result = evidence.verify()
    assert result.verification_status == "passed"
    assert result.science.status == "not_estimable"
    assert result.science.improvement["passed"] is None
    assert all(row["status"] == "not_estimable" for row in result.science.failure_conditions)


@pytest.mark.parametrize("change", ["duplicate_scope", "duplicate_pair", "missing_task", "wrong_gate", "wrong_verifier", "old_root_latest"])
def test_status_and_list_length_do_not_grant_receipt_admission(evidence, change):
    if change == "duplicate_scope":
        evidence.native["native_scopes"][-1] = deepcopy(evidence.native["native_scopes"][-2])
    elif change == "duplicate_pair":
        evidence.native["pairs"][-1] = deepcopy(evidence.native["pairs"][-2])
    elif change == "missing_task":
        evidence.native["native_scopes"][0]["durable_receipts"].pop()
    elif change == "old_root_latest":
        evidence.native["latest_comparisons"][0]["historical_root_comparison_used"] = True
    else:
        if change == "wrong_gate":
            evidence.scientific["gate"]["protocol_sha256"] = "0" * 64
        else:
            evidence.scientific["gate"]["actual_native_verifier_path"] = evidence.source_pin["path"]
        evidence.update_science()
        assert evidence.verify().verification_status == "failed"
        return
    evidence.update_native()
    assert evidence.verify().verification_status == "failed"


@pytest.mark.parametrize("change", ["promoted_hypothesis", "wrong_threshold", "dropped_class", "wrong_denominator", "perfect_degenerate", "omitted_field", "fake_delta"])
def test_scientific_reporting_contract_cannot_be_relaxed(evidence, change):
    primary = evidence.scientific["branches"][1]["reference_fields"][PRIMARY]
    if change == "promoted_hypothesis":
        evidence.scientific["primary_hypothesis"]["status"] = "passed"
    elif change == "wrong_threshold":
        primary["fixed_thresholds"]["adjusted_rand_index_minimum"] = .3
    elif change == "dropped_class":
        primary["coverage"]["per_class"].pop("unknown")
    elif change == "wrong_denominator":
        primary["coverage"]["valid_reference"]["full_source_cells"] = 7
    elif change == "perfect_degenerate":
        evidence.scientific["branches"][0]["reference_fields"][SECONDARY]["agreement"]["adjusted_rand_index"] = 1.0
    elif change == "fake_delta":
        evidence.scientific["primary_hypothesis"]["ari_improvement"] = .2
    else:
        evidence.scientific["branches"][0]["reference_fields"].pop(SECONDARY)
    evidence.update_science()
    assert evidence.verify().verification_status == "failed"


def test_source_change_requires_explicit_bound_archive(evidence):
    archive, sha, archived_source = evidence.archive()
    evidence.source.write_text("P5 evolved current source\n")
    assert evidence.verify().verification_status == "failed"
    result = evidence.verify(archive_receipt_path=archive, expected_archive_receipt_sha256=sha)
    assert result.verification_status == "passed"
    substitution = next(row for row in result.verified_files if row.path == str(evidence.source))
    assert substitution.source == "bound_archive" and substitution.verified_path == str(archived_source)


@pytest.mark.parametrize("change", ["wrong_receipt_sha", "wrong_index_path", "wrong_index_sha", "duplicate_mapping", "missing_mapping", "wrong_source_sha", "changed_archive", "data_replacement"])
def test_archive_cannot_hide_missing_or_unbound_evidence(evidence, change):
    archive, sha, archived_source = evidence.archive()
    value = json.loads(archive.read_text())
    evidence.source.write_text("changed current source\n")
    if change == "wrong_receipt_sha":
        sha = "0" * 64
    elif change == "wrong_index_path":
        value["p4_index"]["path"] = str(evidence.directory / "other.json")
    elif change == "wrong_index_sha":
        value["p4_index"]["sha256"] = "0" * 64
    elif change == "duplicate_mapping":
        value["sources"].append(value["sources"][0])
    elif change == "missing_mapping":
        value["sources"] = []
    elif change == "wrong_source_sha":
        value["sources"][0]["sha256"] = "0" * 64
    elif change == "changed_archive":
        archived_source.write_text("archive bytes changed\n")
    else:
        value["sources"] = [dict(original_path=evidence.data["path"], archive_path=str(archived_source), sha256=evidence.data["sha256"])]
    if change != "wrong_receipt_sha":
        archive.write_text(json.dumps(value))
        sha = digest(archive)
    assert evidence.verify(archive_receipt_path=archive, expected_archive_receipt_sha256=sha).verification_status == "failed"


def test_changed_archive_is_rejected_even_if_current_source_still_matches(evidence):
    archive, sha, archived_source = evidence.archive()
    archived_source.write_text("changed unused copy\n")
    assert evidence.verify(archive_receipt_path=archive, expected_archive_receipt_sha256=sha).verification_status == "failed"


@pytest.mark.parametrize("kind", ["same_size_bytes", "size", "index_sha", "index_duplicate", "pinned_publication_log"])
def test_hash_size_index_and_publication_boundaries(evidence, kind):
    if kind == "same_size_bytes":
        path = Path(evidence.data["path"])
        path.write_bytes(b"x" * path.stat().st_size)
    elif kind == "size":
        evidence.index["evidence_files"][0]["bytes"] += 1
        evidence.index_path.write_text(json.dumps(evidence.index))
        evidence.index_sha = digest(evidence.index_path)
    elif kind == "index_sha":
        evidence.index_sha = "0" * 64
    elif kind == "index_duplicate":
        evidence.index["evidence_files"].append(evidence.index["evidence_files"][0])
        evidence.index_path.write_text(json.dumps(evidence.index))
        evidence.index_sha = digest(evidence.index_path)
    else:
        evidence.write("publication.log", b"published but not closed\n")
        evidence.publish()
    result = evidence.verify()
    assert result.verification_status == "failed"
    assert result.engineering.status == result.science.status == "not_verified"


def test_excluded_log_append_does_not_modify_frozen_closure(evidence):
    with evidence.excluded_log.open("a") as stream:
        stream.write("more publication output\n")
    result = evidence.verify()
    assert result.verification_status == "passed"
    assert str(evidence.excluded_log) in result.publication_logs_outside_index


def test_missing_archive_argument_partner_fails(evidence):
    archive, _, _ = evidence.archive()
    assert evidence.verify(archive_receipt_path=archive).verification_status == "failed"


def test_index_science_copy_must_equal_pinned_checker(evidence):
    evidence.index["scientific_primary_hypothesis"]["status"] = "passed"
    evidence.publish()
    assert evidence.verify().verification_status == "failed"


@pytest.mark.parametrize("change", ["missing_capability", "failed_output", "approximate_array", "nonzero_tolerance", "changed_parameter", "empty_projects"])
def test_regression_top_level_pass_does_not_hide_failed_or_missing_checks(evidence, change):
    value = json.loads(Path(evidence.index["primary_scientific_regression"]["path"]).read_text())
    if change == "empty_projects":
        value["projects"] = []
    elif change == "nonzero_tolerance":
        value["comparison_policy"]["numeric_atol"] = 1e-6
    else:
        capability = value["projects"][0]["capabilities"][0]
        if change == "missing_capability":
            value["projects"][0]["capabilities"].pop()
        elif change == "failed_output":
            capability["outputs"][0]["passed"] = False
        elif change == "approximate_array":
            capability["outputs"][0]["semantic"]["X"]["exact"] = False
        else:
            capability["outputs"][0]["new_effective_parameters"] = {"different": True}
    evidence.index["primary_scientific_regression"] = evidence.write_json("regression.json", value)
    evidence.publish()
    assert evidence.verify().verification_status == "failed"


@pytest.mark.parametrize("kind", ["recovery_copy", "ci_copy", "ci_xml", "current_pointer"])
def test_complete_engineering_reporting_is_bound_to_its_original_receipts(evidence, kind):
    if kind == "recovery_copy":
        evidence.index["fault_recovery"]["pbmc3k"]["runs"] = 6
    elif kind == "ci_copy":
        evidence.index["ci"]["lanes"][0]["tests"] = 99
    elif kind == "ci_xml":
        lane = evidence.index["ci"]["lanes"][0]
        lane["junit"] = evidence.write("ci-core.xml", b'<testsuites><testsuite tests="1" failures="0" errors="0" skipped="0"/></testsuites>')
    else:
        item = evidence.native["latest_comparisons"][0]
        pointer = json.loads(Path(item["pointer_path"]).read_text())
        pointer["comparison"] = evidence.data["path"]
        updated = evidence.write_json("pbmc3k/latest_comparison.json", pointer)
        item["pointer_sha256"] = updated["sha256"]
        evidence.update_native()
        assert evidence.verify().verification_status == "failed"
        return
    evidence.publish()
    assert evidence.verify().verification_status == "failed"


def test_generic_artifact_index_is_explicitly_unsupported(evidence):
    evidence.index = dict(schema_version=1, evidence_files=[])
    evidence.publish()
    result = evidence.verify()
    assert result.verification_status == "failed"
    assert result.findings[-1].code == "unsupported_index_kind"


@pytest.mark.parametrize("change", ["protocol_timing", "supplement_timing", "labels_supplied", "label_masks"])
def test_byte_consistent_but_non_prespecified_freeze_is_rejected(evidence, change):
    supplement_path = next(path for path, sha in evidence.native["protected_files"].items()
                           if sha == evidence.native["supplement_sha256"])
    supplement = json.loads(Path(supplement_path).read_text())
    if change == "protocol_timing":
        evidence.protocol["status"] = "created_after_scores"
        pin = evidence.write_json("protocol.json", evidence.protocol)
        supplement["protocol_sha256"] = evidence.native["protocol_sha256"] = evidence.scientific["gate"]["protocol_sha256"] = pin["sha256"]
    elif change == "label_masks":
        pin = evidence.write_json("masks.json", dict(reference_values_used_for_masks=True))
        supplement["masks_sha256"] = evidence.native["masks_sha256"] = evidence.scientific["gate"]["masks_sha256"] = pin["sha256"]
    else:
        pin = None
        if change == "supplement_timing":
            supplement["status"] = "not_frozen"
        else:
            supplement["reference_labels_supplied"] = True
    if pin:
        evidence.native["protected_files"][pin["path"]] = pin["sha256"]
        evidence.addendum["original_files"][pin["path"]] = pin["sha256"]
    freeze = evidence.write_json("supplement.json", supplement)
    evidence.native["supplement_sha256"] = evidence.scientific["gate"]["supplement_sha256"] = freeze["sha256"]
    evidence.native["protected_files"][freeze["path"]] = freeze["sha256"]
    evidence.addendum["original_files"][freeze["path"]] = freeze["sha256"]
    addendum = evidence.write_json("addendum.json", evidence.addendum)
    evidence.index["post_native_validation_addendum"] = addendum
    evidence.native["validation_addendum_sha256"] = evidence.scientific["gate"]["validation_addendum_sha256"] = addendum["sha256"]
    evidence.update_native()
    result = evidence.verify()
    assert result.verification_status == "failed"
    assert result.findings[-1].code == "freeze_scope"


def test_scoring_adapter_identity_cannot_be_omitted(evidence):
    evidence.scientific["code_identity"]["files"] = {
        evidence.native_code["path"]: evidence.native_code["sha256"]}
    evidence.update_science()
    assert evidence.verify().verification_status == "failed"


@pytest.mark.parametrize("change", ["signature", "artifact_hash", "run_dir", "exposure"])
def test_science_cannot_be_rebound_to_another_native_run_or_exposure(evidence, change):
    if change == "exposure":
        evidence.index["prior_knowledge_and_exposure"] = {"fully_blinded": True}
        evidence.publish()
    else:
        row = evidence.scientific["branches"][0]
        if change == "signature":
            row["native_receipts"][0]["signature"] = "different"
        elif change == "artifact_hash":
            row["native_receipts"][0]["artifact_hashes"] = {"different_uri": "different_sha"}
        else:
            row["run_dir"] = str(evidence.directory / "different_run")
        evidence.update_science()
    assert evidence.verify().verification_status == "failed"
