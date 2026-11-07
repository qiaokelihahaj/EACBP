"""Verify a published acceptance index without running its scientific analysis.

This verifies stored bytes, declared scope and relationships between independent
receipts. Scientific conclusions remain the results of their original checker.
Reference data are hashed, never parsed, and no comparison capability is called.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path, PureWindowsPath
import re
from typing import Any, Iterable
import xml.etree.ElementTree as ET

from eacbp.schemas.verification import (
    ReportedEngineering, ReportedScience, VerificationFinding, VerificationResult,
    VerifiedFile,
)


class _InvalidEvidence(ValueError):
    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def _require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise _InvalidEvidence(code, message)


def _identity(path: str | Path) -> str:
    value = Path(path)
    _require(value.is_absolute(), "file_path", f"Evidence path is not absolute: {path}")
    # Windows paths are case insensitive, including paths on a non-Windows test host.
    text = str(value.resolve())
    return text.casefold() if PureWindowsPath(str(path)).drive else text


def _digest(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _sha(value: Any) -> str:
    _require(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None,
             "file_pin", "Supply a lowercase SHA256 file pin")
    return value


def _json_bytes(data: bytes) -> dict[str, Any]:
    def reject_constant(value: str):
        raise _InvalidEvidence("json_values", f"Nonfinite JSON value: {value}")

    value = json.loads(data.decode("utf-8-sig"), parse_constant=reject_constant)
    _require(isinstance(value, dict), "json_shape", "Evidence JSON must contain an object")
    return value


class _Files:
    """Resolve only indexed files and explicitly bound historical source copies."""

    def __init__(self, index: dict, index_path: Path, index_sha: str,
                 archive_path: str | Path | None, archive_sha: str | None):
        self.entries: dict[str, dict] = {}
        self.verified: dict[str, VerifiedFile] = {}
        self.archives: dict[str, dict] = {}
        self.archive_files: list[dict] = []
        self.archive_pin: dict | None = None
        self.source_keys = {_identity(path) for path in index["source_files"]}
        for item in index["evidence_files"]:
            self._add(item, duplicate_allowed=False)
        for path, digest in index["source_files"].items():
            self._add({"path": path, "sha256": digest}, duplicate_allowed=True)
        excluded = {_identity(path) for path in index.get("publication_logs_outside_index", [])}
        _require(not excluded.intersection(self.entries), "publication_boundary",
                 "A publication log declared outside the index is also indexed")
        _require((archive_path is None) == (archive_sha is None), "archive_binding",
                 "Archive receipt path and explicit SHA must be supplied together")
        if archive_path is not None:
            data = Path(archive_path).read_bytes()
            _require(hashlib.sha256(data).hexdigest() == _sha(archive_sha), "archive_binding",
                     "Archive receipt differs from the explicitly supplied SHA")
            archive = _json_bytes(data)
            binding = archive.get("index", archive.get("p4_index", {}))
            _require(archive.get("status") == "passed"
                     and _identity(binding["path"]) == _identity(index_path)
                     and binding["sha256"] == index_sha, "archive_binding",
                     "Historical source archive is not bound to this acceptance index")
            for item in archive["sources"]:
                key = _identity(item["original_path"])
                _require(key not in self.archives, "archive_binding", "Duplicate archived source identity")
                original_drive = PureWindowsPath(item["original_path"]).drive.casefold()
                index_drive = PureWindowsPath(str(index_path)).drive.casefold()
                eligible = key in self.source_keys or bool(original_drive and original_drive != index_drive)
                _require(key in self.entries and eligible
                         and item["sha256"] == self.entries[key]["sha256"], "archive_binding",
                         "Archive mapping is unindexed, changes the expected SHA, or replaces data evidence")
                _identity(item["archive_path"])
                archived = Path(item["archive_path"])
                _require(archived.is_file() and _digest(archived) == item["sha256"], "archive_binding",
                         "Historical source archive bytes changed")
                self.archive_files.append({"original_path": item["original_path"], "path": str(archived.resolve()),
                                           "sha256": item["sha256"], "bytes": archived.stat().st_size})
                self.archives[key] = item
            self.archive_pin = {"path": str(Path(archive_path).resolve()), "sha256": archive_sha, "bytes": len(data)}

    def _add(self, item: dict, *, duplicate_allowed: bool):
        key = _identity(item["path"])
        _sha(item["sha256"])
        if "bytes" in item:
            _require(type(item["bytes"]) is int and item["bytes"] >= 0,
                     "file_pin", "Indexed file size must be a nonnegative integer")
        if key in self.entries:
            previous = self.entries[key]
            _require(duplicate_allowed and previous["sha256"] == item["sha256"],
                     "file_identity", "Duplicate or conflicting indexed file identity")
        else:
            self.entries[key] = dict(item)

    def verify(self, path: str | Path, digest: str | None = None) -> VerifiedFile:
        key = _identity(path)
        _require(key in self.entries, "file_binding", f"Receipt references an unindexed file: {path}")
        item = self.entries[key]
        _require(digest is None or digest == item["sha256"], "file_binding", f"Receipt and index SHA differ: {path}")
        if key in self.verified:
            return self.verified[key]
        original = Path(item["path"])
        actual = original
        source = "current_file"
        matches = actual.is_file() and _digest(actual) == item["sha256"]
        if not matches and key in self.archives:
            actual = Path(self.archives[key]["archive_path"])
            source = "bound_archive"
            matches = actual.is_file() and _digest(actual) == item["sha256"]
        _require(matches, "file_integrity", f"Missing or changed indexed file without a matching bound archive: {original}")
        size = actual.stat().st_size
        _require("bytes" not in item or size == item["bytes"], "file_integrity", f"Indexed file size differs: {original}")
        result = VerifiedFile(path=str(original), sha256=item["sha256"], bytes=size,
                              verified_path=str(actual.resolve()), source=source)
        self.verified[key] = result
        return result

    def pin(self, value: dict) -> dict:
        result = self.verify(value["path"], _sha(value["sha256"]))
        _require("bytes" not in value or (type(value["bytes"]) is int and value["bytes"] == result.bytes),
                 "file_binding", "Receipt and actual file size differ")
        return {"path": result.path, "sha256": result.sha256, "bytes": result.bytes}

    def read(self, value: dict) -> dict:
        _require(Path(value["path"]).suffix.casefold() == ".json", "json_shape", "Only bound JSON receipts may be parsed")
        result = self.verify(value["path"], _sha(value["sha256"]))
        data = Path(result.verified_path).read_bytes()
        _require(hashlib.sha256(data).hexdigest() == result.sha256 and len(data) == result.bytes,
                 "file_integrity", "Evidence changed between hash verification and JSON reading")
        return _json_bytes(data)

    def by_sha(self, digest: str, candidates: dict | None = None) -> dict:
        _sha(digest)
        matches = [item for key, item in self.entries.items() if item["sha256"] == digest
                   and (candidates is None or any(_identity(path) == key and value == digest
                                                 for path, value in candidates.items()))]
        _require(bool(matches), "file_binding", "Receipt SHA has no indexed file binding")
        return self.pin(matches[0])

    def mapping(self, mapping: dict) -> None:
        _require(isinstance(mapping, dict) and bool(mapping), "code_identity", "Empty file identity map")
        for path, digest in mapping.items():
            self.verify(path, _sha(digest))

    def named(self, name: str) -> dict:
        matches = [item for item in self.entries.values() if Path(item["path"]).name == name]
        _require(len(matches) == 1, "file_binding", f"Expected one indexed {name} receipt")
        return self.pin(matches[0])


def _unique(rows: Iterable[dict], fields: tuple[str, ...], name: str) -> set[tuple]:
    rows = list(rows)
    identities = [tuple(row[field] for field in fields) for row in rows]
    _require(len(identities) == len(set(identities)), "scope_identity", f"Duplicate {name} identity")
    return set(identities)


def _nonnegative(value: Any) -> bool:
    return type(value) is int and value >= 0


def _disposition(row: dict) -> None:
    keys = ("full_source_cells", "branch_raw_cells", "retained_cells",
            "omitted_by_frozen_source_selection", "actual_qc_candidate_cells")
    _require(all(_nonnegative(row[key]) for key in keys), "scientific_report", "Invalid source denominator counts")
    total, raw, kept, omitted, candidates = (row[key] for key in keys)
    _require(raw + omitted == total and kept + candidates == raw, "scientific_report",
             "Source omission and QC exclusion do not partition the reported denominator")
    ratios = {"full_source_retention": kept / total if total else None,
              "branch_raw_coverage": raw / total if total else None,
              "qc_retention_among_branch_raw": kept / raw if raw else None}
    _require(all(row[key] == value for key, value in ratios.items()), "scientific_report",
             "Reported retention does not match its stated denominator")


def _field_result(field: dict, policy: dict) -> None:
    c, a = field["coverage"], field["agreement"]
    for name in ("raw_source_n", "valid_label_n", "missing_label_n"):
        _require(_nonnegative(c[name]), "scientific_report", "Invalid source or label count")
    _require(c["valid_label_n"] + c["missing_label_n"] == c["raw_source_n"],
             "scientific_report", "Valid and missing labels do not cover the full source")
    for key, total in (("full_source", c["raw_source_n"]), ("valid_reference", c["valid_label_n"]),
                       ("missing_reference", c["missing_label_n"])):
        _disposition(c[key])
        _require(c[key]["full_source_cells"] == total, "scientific_report", "Source denominator changed")
    for row in c["per_class"].values():
        _disposition(row)
        _require(row["full_source_cells"] > 0, "scientific_report", "An empty class was invented")
    for key in ("full_source_cells", "branch_raw_cells", "retained_cells", "omitted_by_frozen_source_selection", "actual_qc_candidate_cells"):
        _require(sum(row[key] for row in c["per_class"].values()) == c["valid_reference"][key],
                 "scientific_report", "Per-class denominators do not cover the valid-label source")
        _require(c["valid_reference"][key] + c["missing_reference"][key] == c["full_source"][key],
                 "scientific_report", "Valid and missing dispositions do not cover the full source")
    _require(c["all_valid_classes_kept_in_denominator"] is True, "scientific_report", "Reference classes omitted")
    for key in ("common_valid_label_cells", "reference_common_classes", "prediction_common_classes"):
        _require(_nonnegative(a[key]), "scientific_report", "Invalid common-cell/class count")
    _require(a["common_valid_label_cells"] == c["valid_reference"]["retained_cells"],
             "scientific_report", "Metric population differs from retained valid-label cells")
    _require(a["reference_common_classes"] == sum(row["retained_cells"] > 0 for row in c["per_class"].values())
             and a["prediction_common_classes"] <= a["common_valid_label_cells"],
             "scientific_report", "Reported common classes contradict the retained reference denominator")
    estimable = a["common_valid_label_cells"] >= 2 and a["reference_common_classes"] >= 2 and a["prediction_common_classes"] >= 2
    _require(a["status"] == ("estimable" if estimable else "not_estimable") and a["nmi_average_method"] == "arithmetic",
             "scientific_report", "Reported estimability or NMI convention differs")
    ari, nmi = a["adjusted_rand_index"], a["normalized_mutual_information"]
    _require((estimable and all(type(x) in (int, float) and math.isfinite(x) for x in (ari, nmi))
              and -1 <= ari <= 1 and 0 <= nmi <= 1) or (not estimable and ari is None and nmi is None),
             "scientific_report", "Degenerate metrics must be null; estimable metrics must be finite")
    thresholds = {"adjusted_rand_index_minimum": policy["minimum_absolute_ari"],
                  "normalized_mutual_information_minimum": policy["minimum_nmi"],
                  "retention_minimum": policy["minimum_overall_labelled_source_retention"]}
    _require(field["fixed_thresholds"] == thresholds, "scientific_report", "Reported frozen thresholds changed")
    overall = c["valid_reference"]["full_source_retention"]
    expected = {"overall_valid_reference_retention": overall is not None and overall >= policy["minimum_overall_labelled_source_retention"],
                "every_valid_reference_class_retention": bool(c["per_class"]) and all(row["full_source_retention"] >= policy["minimum_per_class_source_retention"] for row in c["per_class"].values()),
                "adjusted_rand_index": None if ari is None else ari >= policy["minimum_absolute_ari"],
                "normalized_mutual_information": None if nmi is None else nmi >= policy["minimum_nmi"]}
    _unique(field["checks"], ("name",), "scientific check")
    actual = {row["name"]: row["passed"] for row in field["checks"]}
    _require(actual == expected and all(value is None or type(value) is bool for value in actual.values()),
             "scientific_report", "Scientific check does not follow the frozen reported metric/denominator")
    _require(field["failed_checks"] == [row["name"] for row in field["checks"] if row["passed"] is False]
             and field["not_estimable_checks"] == [row["name"] for row in field["checks"] if row["passed"] is None],
             "scientific_report", "Failure or not-estimable reasons were omitted")
    status = "not_estimable" if not estimable else "passed" if all(expected.values()) else "failed"
    _require(field["status"] == status, "scientific_report", "Scientific status contradicts its reported checks")


def _regression_checks(value: Any) -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            if key in {"passed", "exact", "keys_unchanged"}:
                _require(item is True, "engineering_receipt", "Regression contains a failed or approximate scientific check")
            if key in {"rtol", "atol"}:
                _require(type(item) in (int, float) and item == 0, "engineering_receipt", "Regression numeric tolerance changed")
            _regression_checks(item)
    elif isinstance(value, list):
        for item in value:
            _regression_checks(item)


def _verify_reports(index: dict, files: _Files) -> tuple[dict, dict, dict, dict, dict]:
    native_pin = files.pin(index["independent_native_gate"])
    science_pin = files.pin(index["external_science_gate"])
    addendum_pin = files.pin(index["post_native_validation_addendum"])
    failed_pin = files.pin(index["preserved_failed_native_gate"])
    native, scientific = files.read(native_pin), files.read(science_pin)
    addendum, failed = files.read(addendum_pin), files.read(failed_pin)
    _require(native["status"] == "passed" and native["all_10_native_scopes_admitted_and_question_complete"] is True
             and native["all_8_pairs_independently_recomputed"] is True
             and native["reference_labels_opened"] is False and native["external_scores_opened"] is False,
             "engineering_receipt", "Native independent checker did not close its declared scope")
    protected = native["protected_files"]
    supplement_pin = files.by_sha(native["supplement_sha256"], protected)
    protocol_pin = files.by_sha(native["protocol_sha256"], protected)
    masks_pin = files.by_sha(native["masks_sha256"], protected)
    supplement, protocol, masks = files.read(supplement_pin), files.read(protocol_pin), files.read(masks_pin)
    _require(supplement["status"] == "frozen_before_first_native_execution"
             and supplement["reference_labels_supplied"] is False
             and protocol["status"] == "prespecified_before_native_and_external_scores"
             and masks["reference_values_used_for_masks"] is False,
             "freeze_scope", "Freeze/protocol/masks do not preserve the declared reference-free timing")
    _require(supplement["protocol_sha256"] == protocol_pin["sha256"] and supplement["masks_sha256"] == masks_pin["sha256"],
             "gate_binding", "Native scope uses another freeze/protocol/mask")
    files.pin({"path": supplement["protocol_path"], "sha256": supplement["protocol_sha256"]})
    files.pin({"path": supplement["masks_path"], "sha256": supplement["masks_sha256"]})
    expected_scopes = set()
    projects = set()
    for dataset in supplement["datasets"]:
        project = dataset["project"]
        _require(project not in projects, "scope_identity", "Duplicate frozen project")
        projects.add(project)
        prepared = dataset["prepared"]
        bindings = files.read({"path": prepared["bindings_path"], "sha256": prepared["bindings_sha256"]})
        ids = _unique(bindings["branches"], ("branch_id",), "frozen branch")
        _require(len(ids) == dataset["expected_branch_count"], "scope_identity", "Frozen branch count differs")
        expected_scopes.update((project, value[0]) for value in ids)
        files.pin({"path": dataset["original_source"], "sha256": dataset["original_source_sha256"]})
    actual_scopes = _unique(native["native_scopes"], ("project", "branch_id"), "native scope")
    baseline = protocol["new_external_check"]["baseline_branch"]
    actual_pairs = _unique(native["pairs"], ("project", "branch_id"), "independent pair")
    _require(len(expected_scopes) == 10 and len(actual_pairs) == 8
             and actual_scopes == expected_scopes and actual_pairs == {row for row in expected_scopes if row[1] != baseline},
             "scope_identity", "Native/pair identities do not cover the exact frozen scope")
    _require(index["native_scopes"] == native["native_scopes"] and index["independent_pairs"] == native["pairs"],
             "report_binding", "Index changed independently reported native/pair results")
    for scope in native["native_scopes"]:
        _require(scope["native_status"] == "success" and scope["question_status"] == "complete"
                 and scope["native_payload_values_loaded_only_through_get_audited"] is True,
                 "engineering_receipt", "Native scope lacks reported durable admission")
        _unique(scope["durable_receipts"], ("task_id",), "durable task")
        _require({row["capability"] for row in scope["durable_receipts"]} == {
            "dataset_audit", "qc_review", "qc", "normalization", "integration", "clustering", "cluster_markers"}
            and len(scope["durable_receipts"]) == 7
            and all(row["get_audited_verified"] is True for row in scope["durable_receipts"]),
                 "engineering_receipt", "Native task evidence is missing or not admitted")
        for task in scope["durable_receipts"]:
            _sha(task["signature"])
            _require(bool(task["artifact_hashes"]), "engineering_receipt", "Admitted task has no output artifact hashes")
            for artifact_digest in task["artifact_hashes"].values():
                _require(isinstance(artifact_digest, str) and artifact_digest.startswith("sha256:"),
                         "file_binding", "Artifact receipt must declare its serialized SHA256")
                files.by_sha(artifact_digest.removeprefix("sha256:"))
    for receipt in (native, scientific):
        files.mapping(receipt["protected_files"])
        files.mapping(receipt["code_identity"]["files"])
    files.mapping(supplement["checker_files"])
    _require(protocol["checker_files"] == supplement["checker_files"], "code_identity", "Frozen checker maps differ")
    _require(addendum["status"] == "validation_only_addendum_post_native_pre_external_scores"
             and addendum["science_hypothesis_and_native_configuration_unchanged"] is True
             and type(addendum["native_rerun_count"]) is int and addendum["native_rerun_count"] == 0
             and addendum["sdrf_reference_values_read_by_score_evaluator"] is False,
             "validation_history", "Validation addendum changed its declared scientific/exposure scope")
    files.mapping(addendum["original_files"])
    files.mapping(addendum["new_files"])
    original_keys, new_keys = ({_identity(path) for path in addendum[name]} for name in ("original_files", "new_files"))
    required_original = {_identity(pin["path"]) for pin in (supplement_pin, protocol_pin, masks_pin, failed_pin)}
    required_original.update(_identity(path) for path in supplement["checker_files"])
    _require(required_original <= original_keys and not original_keys.intersection(new_keys),
             "validation_history", "Original freeze/checker/failure pins missing or relabelled as new files")
    _require(files.pin(addendum["original_failed_gate"]) == failed_pin
             and failed["status"] == "failed_native_scope" and failed["reference_labels_opened"] is False
             and failed["external_scores_opened"] is False,
             "validation_history", "Original failed gate or its exposure state was rewritten")
    files.pin(addendum["diagnosis"])
    _require(_identity(addendum["diagnosis"]["path"]) in original_keys, "validation_history", "Unbound numerical diagnosis")
    _require(bool(addendum["verification_receipts"]), "validation_history", "Addendum has no verification receipts")
    for pin in addendum["verification_receipts"]:
        files.pin(pin)
    gate = scientific["gate"]
    _require(files.pin(gate) == native_pin
             and gate["validation_addendum_sha256"] == native["validation_addendum_sha256"] == addendum_pin["sha256"]
             and _identity(gate["validation_addendum_path"]) == _identity(addendum_pin["path"])
             and _identity(native["validation_addendum_path"]) == _identity(addendum_pin["path"])
             and gate["protocol_sha256"] == protocol_pin["sha256"] and gate["supplement_sha256"] == supplement_pin["sha256"]
             and gate["masks_sha256"] == masks_pin["sha256"] and gate["environment"] == supplement["environment"],
             "gate_binding", "Scientific receipt is not bound to this native gate/addendum/freeze")
    verifier_key = _identity(gate["actual_native_verifier_path"])
    _require(verifier_key in new_keys and verifier_key in {_identity(path) for path in native["code_identity"]["files"]},
             "code_identity", "Actual native verifier is absent from its new-file/code identity binding")
    scientific_new_keys = {_identity(path) for path in scientific["code_identity"]["files"]}.intersection(new_keys) - {verifier_key}
    _require(bool(scientific_new_keys), "code_identity", "Scientific code identity lacks a new versioned scoring adapter")
    _require(gate["all_10_native_scopes_currently_get_audited"] is True and gate["all_8_pairs_independently_recomputed"] is True
             and scientific["status"] == "evaluated" and scientific["engineering_status"] == "passed"
             and scientific["reference_labels_opened_after_all_native_scopes_and_pairs"] is True,
             "engineering_receipt", "External checker lacks a completed native-before-reference gate")
    _require(index["external_science"] == scientific["branches"]
             and index["scientific_primary_hypothesis"] == scientific["primary_hypothesis"]
             and index["external_reference"] == scientific["reference"]
             and index["prior_knowledge_and_exposure"] == scientific["prior_knowledge_and_exposure"]
             and index["scientific_interpretation"] == scientific["interpretation"]
             and index["biological_validation_established"] is False,
             "report_binding", "Index changed original scientific results/reference/exposure or promoted biological validation")
    for receipt in (index, scientific):
        _require(all(receipt[name] is False for name in ("winner_selected", "default_changed", "primary_replaced")),
                 "scientific_report", "Verification cannot select a winner or replace the primary")
    policy = protocol["new_external_check"]
    fields = {policy["primary_field"], policy["secondary_field"]}
    expected_science = {row for row in expected_scopes if row[0] == policy["project"]}
    _require(_unique(scientific["branches"], ("project", "branch_id"), "external branch") == expected_science,
             "scope_identity", "External checker did not report every frozen branch")
    for branch in scientific["branches"]:
        native_scope = next(row for row in native["native_scopes"]
                            if (row["project"], row["branch_id"]) == (branch["project"], branch["branch_id"]))
        _unique(branch["native_receipts"], ("task_id",), "scientific native task")
        identity_keys = ("task_id", "capability", "signature", "artifact_hashes")
        scientific_tasks = {row["task_id"]: {key: row[key] for key in identity_keys} for row in branch["native_receipts"]}
        native_tasks = {row["task_id"]: {key: row[key] for key in identity_keys} for row in native_scope["durable_receipts"]}
        _require(scientific_tasks == native_tasks
                 and all(row["all_siblings_get_audited"] is True for row in branch["native_receipts"])
                 and _identity(branch["run_dir"]) == _identity(native_scope["run_dir"]),
                 "gate_binding", "Scientific branch uses different native task/artifact identities or run directory")
        _require(set(branch["reference_fields"]) == fields, "scientific_report", "Preselected reference fields omitted or substituted")
        for field in branch["reference_fields"].values():
            _field_result(field, policy)
    reference = scientific["reference"]
    _require(reference["canonical_sha256"] == policy["reference"]["sha256"]
             and _identity(reference["canonical_path"]) == _identity(policy["reference"]["path"])
             and reference["fields"] == {"primary": policy["primary_field"], "secondary_descriptive": policy["secondary_field"]}
             and set(reference["missing_tokens_stripped_casefold"]) == set(policy["missing_tokens_casefold"]),
             "reference_binding", "Reported reference/field/missing-token source changed")
    files.pin(policy["reference"])
    rows = {row["branch_id"]: row for row in scientific["branches"]}
    target = rows[policy["hypothesis_branch"]]["reference_fields"][policy["primary_field"]]
    hypothesis = scientific["primary_hypothesis"]
    improvement = hypothesis["improvement_check"]
    delta = hypothesis["ari_improvement"]
    first = rows[baseline]["reference_fields"][policy["primary_field"]]["agreement"]["adjusted_rand_index"]
    second = target["agreement"]["adjusted_rand_index"]
    expected_delta = second - first if first is not None and second is not None else None
    expected_improvement = "not_estimable" if delta is None else "passed" if delta >= policy["minimum_ari_improvement"] else "failed"
    expected_combined = "not_estimable" if expected_improvement == "not_estimable" or target["status"] == "not_estimable" else "passed" if expected_improvement == target["status"] == "passed" else "failed"
    _require(hypothesis["reference_field"] == policy["primary_field"]
             and hypothesis["primary_ari"] == first and hypothesis["resolution_0_5_ari"] == second and delta == expected_delta
             and hypothesis["minimum_improvement"] == policy["minimum_ari_improvement"]
             and improvement["status"] == expected_improvement and improvement["ari_improvement"] == delta
             and improvement["minimum_improvement"] == policy["minimum_ari_improvement"]
             and improvement["passed"] is (None if delta is None else delta >= policy["minimum_ari_improvement"])
             and hypothesis["status"] == expected_combined
             and hypothesis["target_acceptance"]["branch_id"] == policy["hypothesis_branch"]
             and hypothesis["target_acceptance"]["reference_field"] == policy["primary_field"]
             and hypothesis["target_acceptance"]["status"] == target["status"]
             and hypothesis["target_acceptance"]["checks"] == target["checks"]
             and hypothesis["target_acceptance"]["failed_checks"] == target["failed_checks"]
             and hypothesis["target_acceptance"]["not_estimable_checks"] == target["not_estimable_checks"]
             and hypothesis["target_acceptance"]["thresholds"] == target["fixed_thresholds"],
             "scientific_report", "Combined hypothesis was promoted by improvement alone or changed its scope")
    regression = files.read(index["primary_scientific_regression"])
    _require(regression["status"] == "passed" and regression["both_pbmc_primary_native_scopes_and_durable_audits_verified"] is True,
             "engineering_receipt", "Primary scientific-array regression did not pass its admission gate")
    capabilities = {"qc", "normalization", "integration", "clustering", "cluster_markers"}
    comparison_policy = regression["comparison_policy"]
    _require(set(comparison_policy["required_capabilities"]) == capabilities
             and all(comparison_policy[key] is True for key in ("raw_count_layers_and_dtype_exact", "scientific_arrays_and_tables_exact", "cell_and_feature_order_exact", "method_operation_effective_scientific_parameters_must_match"))
             and all(comparison_policy[key] == 0 for key in ("numeric_rtol", "numeric_atol", "embedding_rtol", "embedding_atol")),
             "engineering_receipt", "Primary regression changed its zero-tolerance scientific policy")
    _require(_unique(regression["projects"], ("project",), "regression project") == {("pbmc3k",)},
             "scope_identity", "Primary regression is not the declared PBMC-only scope")
    project = regression["projects"][0]
    _require(project["status"] == "passed" and project["failed_capabilities"] == []
             and _unique(project["capabilities"], ("capability",), "regression capability") == {(name,) for name in capabilities},
             "engineering_receipt", "Primary regression omitted or failed a capability")
    for capability in project["capabilities"]:
        _require(capability["passed"] is True and capability["old_method"] == capability["new_method"]
                 and bool(capability["outputs"]), "engineering_receipt", "Primary regression method/output check failed")
        for output in capability["outputs"]:
            _require(output["passed"] is True and output["effective_scientific_parameters_equal"] is True
                     and output["old_effective_parameters"] == output["new_effective_parameters"]
                     and output["semantic"]["passed"] is True,
                     "engineering_receipt", "Primary regression output changed its reported scientific values")
            _regression_checks(output["semantic"])
    files.mapping(regression["code_identity"]["source_files"])
    for name in ("protected_old_files", "protected_new_files"):
        files.mapping(regression[name])
    recovery_pin = files.pin(index["recovery_gate"]) if "recovery_gate" in index else files.named("p4-independent-recovery-v1.json")
    recovery = files.read(recovery_pin)
    _require(recovery["status"] == "passed_independent_recovery"
             and recovery["all_expected_fault_and_recovery_phases_verified"] is True
             and recovery["native_execution_attempts_total"] == len(expected_scopes)
             and recovery["references_opened"] is False and set(recovery["projects"]) == projects
             and recovery["projects"] == index["fault_recovery"]
             and files.pin(recovery["supplement"]) == supplement_pin,
             "engineering_receipt", "Recovery summary differs from its independent frozen-scope receipt")
    files.pin(recovery["checker"])
    _require(supplement["checker_files"].get(recovery["checker"]["path"]) == recovery["checker"]["sha256"],
             "code_identity", "Recovery checker differs from the frozen checker universe")
    lane_ids = _unique(index["ci"]["lanes"], ("lane",), "CI lane")
    _require(lane_ids == {("core",), ("standard",), ("advanced-statistics",)}, "engineering_receipt", "CI lane omitted or duplicated")
    for lane in index["ci"]["lanes"]:
        receipt = files.read(lane)
        _require(receipt["returncode"] == 0 and receipt["failures"] == receipt["errors"] == receipt["skipped"] == 0
                 and receipt["tests"] == lane["tests"] and receipt["explicit_deselect"] == lane["explicit_deselect"]
                 and receipt["environment_before"] == receipt["environment_after"] == supplement["environment"],
                 "engineering_receipt", "CI summary changed its run, exclusions or environment")
        junit = files.verify(lane["junit"]["path"], lane["junit"]["sha256"])
        files.pin(lane["junit"])
        files.pin(lane["log"])
        _require(_identity(receipt["junit_path"]) == _identity(lane["junit"]["path"])
                 and _identity(receipt["log_path"]) == _identity(lane["log"]["path"]), "file_binding", "CI XML/log belongs to another run")
        xml = Path(junit.verified_path).read_bytes()
        _require(hashlib.sha256(xml).hexdigest() == junit.sha256, "file_integrity", "CI XML changed during reading")
        suites = list(ET.fromstring(xml).iter("testsuite"))
        _require(bool(suites) and all(sum(int(suite.attrib[key]) for suite in suites) == receipt[key]
                                     for key in ("tests", "failures", "errors", "skipped")),
                 "engineering_receipt", "CI XML totals do not match the recorded receipt")
    _require(index["ci"]["total_passed"] == sum(lane["tests"] for lane in index["ci"]["lanes"]),
             "engineering_receipt", "Total CI count differs from its distinct lanes")
    latest = native.get("latest_comparisons", [])
    _require(_unique(latest, ("project",), "latest comparison") == {(project,) for project in projects},
             "latest_binding", "Latest comparison delivery omitted a project")
    for item in latest:
        pointer = files.read({"path": item["pointer_path"], "sha256": item["pointer_sha256"]})
        for kind in ("comparison", "report"):
            files.pin({"path": item[kind + "_path"], "sha256": item[kind + "_sha256"]})
            selected = pointer.get(kind + "_path", pointer.get(kind))
            _require(selected is not None and _identity(selected) == _identity(item[kind + "_path"])
                     and pointer[kind + "_sha256"] == item[kind + "_sha256"],
                     "latest_binding", "Latest pointer selects a different historical result")
        _require(pointer["invocation_id"] == item["invocation_id"] and item["historical_root_comparison_used"] is False,
                 "latest_binding", "Historical root comparison substituted for current delivery")
    evidence = {"native_gate": native_pin, "scientific_checker": science_pin, "validation_addendum": addendum_pin,
                "preserved_failed_gate": failed_pin, "supplement": supplement_pin, "protocol": protocol_pin,
                "masks": masks_pin, "reference": files.pin(policy["reference"]), "recovery_checker": recovery_pin,
                "native_code": native["code_identity"], "scientific_code": scientific["code_identity"],
                "original_sources": [{"path": row["original_source"], "sha256": row["original_source_sha256"]} for row in supplement["datasets"]]}
    return native, scientific, policy, evidence, {"preserved_failed_gate": failed_pin, "validation_addendum": addendum_pin,
        "current_native_gate": native_pin, "current_scientific_checker": science_pin, "latest_comparisons": latest}


def _science_summary(scientific: dict, policy: dict) -> ReportedScience:
    hypothesis = scientific["primary_hypothesis"]
    target = next(row for row in scientific["branches"] if row["branch_id"] == policy["hypothesis_branch"])["reference_fields"][policy["primary_field"]]
    observed = {"adjusted_rand_index": target["agreement"]["adjusted_rand_index"],
                "normalized_mutual_information": target["agreement"]["normalized_mutual_information"],
                "overall_valid_reference_retention": target["coverage"]["valid_reference"]["full_source_retention"],
                "every_valid_reference_class_retention": {name: value["full_source_retention"] for name, value in target["coverage"]["per_class"].items()}}
    thresholds = {"adjusted_rand_index": policy["minimum_absolute_ari"], "normalized_mutual_information": policy["minimum_nmi"],
                  "overall_valid_reference_retention": policy["minimum_overall_labelled_source_retention"],
                  "every_valid_reference_class_retention": policy["minimum_per_class_source_retention"]}
    reasons = [{"condition": row["name"], "status": "not_estimable" if row["passed"] is None else "failed",
                "observed": observed[row["name"]], "minimum": thresholds[row["name"]],
                "branch_id": policy["hypothesis_branch"], "reference_field": policy["primary_field"]}
               for row in target["checks"] if row["passed"] is not True]
    if hypothesis["improvement_check"]["status"] != "passed":
        reasons.insert(0, {"condition": "ari_improvement", "status": hypothesis["improvement_check"]["status"],
                           "observed": hypothesis["ari_improvement"], "minimum": policy["minimum_ari_improvement"]})
    return ReportedScience(status=hypothesis["status"], reported_status=hypothesis["status"],
        reference_field=policy["primary_field"], hypothesis=hypothesis, improvement=hypothesis["improvement_check"],
        absolute_acceptance=hypothesis["target_acceptance"], failure_conditions=reasons,
        branches=scientific["branches"], reference=scientific["reference"], interpretation=scientific["interpretation"])


def verify_acceptance_index(index_path: str | Path, *, expected_index_sha256: str,
                            archive_receipt_path: str | Path | None = None,
                            expected_archive_receipt_sha256: str | None = None) -> VerificationResult:
    """Return byte/binding verification and original conclusions; never write files.

    ``expected_index_sha256`` is supplied by the caller. Source evolution requires
    an equally explicit archive receipt pin bound to this exact index. Scientific
    failure or non-estimability is a valid reported outcome, not a closure failure.
    """
    path = Path(index_path).resolve()
    result = VerificationResult(index={"path": str(path), "sha256": expected_index_sha256})
    files = None
    try:
        _sha(expected_index_sha256)
        data = path.read_bytes()
        _require(hashlib.sha256(data).hexdigest() == expected_index_sha256, "index_pin", "Acceptance index SHA differs from the caller's pin")
        index = _json_bytes(data)
        result.index["bytes"] = len(data)
        result.engineering.reported_status = index.get("engineering_status")
        result.science.reported_status = index.get("scientific_primary_hypothesis", {}).get("status")
        _require(index.get("schema_version") == 1, "index_schema", "Unsupported acceptance index schema")
        _require({"independent_native_gate", "external_science_gate", "scientific_primary_hypothesis", "native_scopes", "external_science"} <= set(index),
                 "unsupported_index_kind", "This entry verifies the native/scientific acceptance envelope; a generic artifact index needs its own declared contract")
        files = _Files(index, path, expected_index_sha256, archive_receipt_path, expected_archive_receipt_sha256)
        for item in files.entries.values():
            files.verify(item["path"])
        result.findings.append(VerificationFinding(code="file_closure", status="passed",
            message=f"Verified {len(files.verified)} indexed/source file identities; archive substitutions remain explicit."))
        native, scientific, policy, evidence, relationships = _verify_reports(index, files)
        _require(index["engineering_status"] == "passed", "engineering_receipt", "Index does not report a completed engineering acceptance")
        _require(_digest(path) == expected_index_sha256, "index_pin", "Acceptance index changed during verification")
        if files.archive_pin:
            _require(_digest(files.archive_pin["path"]) == files.archive_pin["sha256"], "archive_binding", "Archive receipt changed during verification")
            evidence["historical_source_archive"] = files.archive_pin
            evidence["archived_source_files"] = files.archive_files
        result.engineering = ReportedEngineering(status="passed", reported_status="passed",
            native_scope_count=len(native["native_scopes"]), independent_pair_count=len(native["pairs"]),
            durable_task_count=sum(len(row["durable_receipts"]) for row in native["native_scopes"]))
        result.science = _science_summary(scientific, policy)
        result.evidence, result.relationships = evidence, relationships
        result.prior_knowledge_and_exposure = index["prior_knowledge_and_exposure"]
        result.publication_logs_outside_index = index.get("publication_logs_outside_index", [])
        result.findings.append(VerificationFinding(code="receipt_closure", status="passed",
            message="Exact frozen scopes, source/code/gate bindings and reported scientific reasons agree. Scientific scores were not re-established."))
        result.verification_status = "passed"
    except (OSError, ValueError, KeyError, TypeError, IndexError, ET.ParseError) as error:
        result.findings.append(VerificationFinding(code=getattr(error, "code", "evidence_structure"),
            status="failed", message=str(error)))
    if files:
        result.verified_files = list(files.verified.values())
    result.summary = format_verification_summary(result)
    return result


def format_verification_summary(result: VerificationResult) -> str:
    """Compact conclusion, reason and pinned evidence for people and CLI output."""
    lines = [f"Evidence verification: {result.verification_status}.",
             f"Engineering: {result.engineering.status}; reported scientific hypothesis: {result.science.status}."]
    if result.verification_status == "failed":
        lines.extend(f"Reason: {item.message}" for item in result.findings if item.status == "failed")
    else:
        lines.append(f"Improvement: {result.science.improvement['status']}; absolute acceptance: {result.science.absolute_acceptance['status']}.")
        target_id = result.science.absolute_acceptance["branch_id"]
        target = next(row for row in result.science.branches if row["branch_id"] == target_id)
        coverage = target["reference_fields"][result.science.reference_field]["coverage"]
        lines.append(f"Reference: {result.science.reference_field}; full source={coverage['raw_source_n']}, valid labels={coverage['valid_label_n']}, missing labels={coverage['missing_label_n']}, classes={len(coverage['per_class'])}.")
        for reason in result.science.failure_conditions:
            lines.append(f"Reason: {reason['condition']} {reason['status']}; observed={reason['observed']}, minimum={reason['minimum']}.")
        for field, value in target["reference_fields"].items():
            if field != result.science.reference_field and value["status"] == "not_estimable":
                lines.append(f"Secondary reference: {field}; not_estimable, ARI/NMI=null, common reference classes={value['agreement']['reference_common_classes']}.")
        archived = sum(row.source == "bound_archive" for row in result.verified_files)
        if archived:
            lines.append(f"Historical source verification used {archived} explicitly bound archive copies; current source may have evolved.")
        lines.append("Scientific verdict is retained from the pinned checker; reference values were not parsed and scores were not recomputed.")
    lines.append(f"Evidence index: {result.index['path']} (SHA256 {result.index['sha256']}).")
    if "scientific_checker" in result.evidence:
        pin = result.evidence["scientific_checker"]
        lines.append(f"Scientific receipt: {pin['path']} (SHA256 {pin['sha256']}).")
    return "\n".join(lines)
