"""Independent acceptance must detect attrition and reversed biological direction."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scripts.evaluate_classic_autonomy import (
    _completed_native_snapshot, _execution_boundary, _independent_donor_effects, _source_subset, clustering_comparison,
    counts_equal, evaluate, leaveout_integrity, population_retention, velocity_direction,
    sha256,
)


def test_rare_class_loss_fails_despite_perfect_remaining_clusters():
    truth = pd.Series(["major"] * 98 + ["rare"] * 2, index=[f"c{i}" for i in range(100)])
    prediction = pd.Series(["0"] * 98, index=truth.index[:98])
    thresholds = {"overall_reference_retention_min": .95, "per_class_retention_min": .8,
                  "adjusted_rand_index_min": .4, "normalized_mutual_information_min": .65}
    result = clustering_comparison(truth, prediction, thresholds)
    assert result["adjusted_rand_index"] == 1
    assert result["overall_retention"] == .98
    assert result["per_class"]["rare"]["retention"] == 0
    assert not next(row for row in result["checks"] if row["name"] == "every_reference_class_retained")["passed"]


def test_duplicate_reference_coverage_is_rejected():
    truth = pd.Series(["a", "b"], index=["one", "two"])
    with pytest.raises(ValueError, match="Repeated"):
        population_retention(truth, ["one", "one"], {})


def test_references_are_not_opened_for_a_running_or_partial_native_analysis(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    # A nonexistent reference directory would fail first if answers were opened
    # before checking the completed native execution boundary.
    for payload in ({"status": "running"}, {"status": "success", "analysis_request_status": "partial"}):
        (run / "summary.json").write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(ValueError, match="must complete before answer-bearing"):
            evaluate("pbmc3k", run, tmp_path / "references-not-present")


def test_velocity_graph_direction_detects_reversal_without_model_rerun():
    labels = pd.Series(["early", "early", "late", "late"])
    graph = sparse.csr_matrix(([2., 3.], ([0, 1], [2, 3])), shape=(4, 4))
    thresholds = {"minimum_cross_state_edges": 2, "direction_index_min": .1}
    forward = velocity_direction(graph, labels, [["early", "late"]], thresholds)
    reversed_result = velocity_direction(graph.T, labels, [["early", "late"]], thresholds)
    assert forward["edges"][0]["direction_index"] == 1
    assert forward["checks"][0]["passed"]
    assert reversed_result["edges"][0]["direction_index"] == -1
    assert not reversed_result["checks"][0]["passed"]


def test_velocity_no_cross_state_evidence_does_not_pass():
    labels = pd.Series(["early", "early", "late", "late"])
    result = velocity_direction(sparse.eye(4), labels, [["early", "late"]],
                                {"minimum_cross_state_edges": 1, "direction_index_min": .1})
    assert not result["checks"][0]["passed"]


@pytest.mark.parametrize("duplicate", [False, True])
def test_stored_zero_or_duplicate_entries_do_not_invent_cross_state_edges(duplicate):
    # Construct noncanonical CSR directly; COO conversion would merge duplicates.
    graph = sparse.csr_matrix((np.array([.5, .5] if duplicate else [1., 0.]),
                               np.array([2, 2] if duplicate else [2, 3]),
                               np.array([0, 2, 2, 2, 2])), shape=(4, 4))
    result = velocity_direction(graph, pd.Series(["early", "early", "late", "late"]),
                                [["early", "late"]], {"minimum_cross_state_edges": 2, "direction_index_min": .1})
    assert result["edges"][0]["direction_index"] == 1
    assert result["edges"][0]["cross_state_edges"] == 1
    assert not result["checks"][0]["passed"]
    assert graph.nnz == 2  # Canonicalization must not alter the caller's graph.


def _admitted_snapshot(summary=None):
    return SimpleNamespace(
        summary=summary or {"status": "success", "analysis_request_status": "complete"},
        task_results=[SimpleNamespace(task_id="native", status=SimpleNamespace(value="success"))],
        audit_contexts={"native": SimpleNamespace(task_id="native", receipt_present=True, status="passed",
                                                   overall_passed=True, stop_rule_triggered=False)},
    )


def test_equivalent_summary_timezone_renderings_are_accepted():
    summary = {"status": "success", "analysis_request_status": "complete", "time": "2026-10-01T07:00:00+00:00"}
    snapshot = _admitted_snapshot({**summary, "time": "2026-10-01T15:00:00+08:00"})
    _completed_native_snapshot(summary, snapshot)


@pytest.mark.parametrize("fault", ["summary", "task", "receipt", "audit", "stop"])
def test_fabricated_summary_or_missing_native_admission_cannot_open_references(tmp_path, monkeypatch, fault):
    import eacbp.artifact.registry
    import eacbp.evidence.snapshot
    run = tmp_path / "run"
    (run / "artifacts").mkdir(parents=True)
    (run / "artifacts" / ".artifact_registry.json").write_text("{}", encoding="utf-8")
    summary = {"status": "success", "analysis_request_status": "complete"}
    (run / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    snapshot = _admitted_snapshot()
    if fault == "summary":
        snapshot.summary = {"status": "failed", "analysis_request_status": "partial"}
    elif fault == "task":
        snapshot.task_results[0].status.value = "failed"
    elif fault == "receipt":
        snapshot.audit_contexts["native"].receipt_present = False
    elif fault == "audit":
        snapshot.audit_contexts["native"].overall_passed = False
    else:
        snapshot.audit_contexts["native"].stop_rule_triggered = True
    order = []
    monkeypatch.setattr(eacbp.evidence.snapshot, "load_study_snapshot", lambda path: snapshot)
    monkeypatch.setattr(eacbp.artifact.registry, "ArtifactRegistry", lambda path: object())
    monkeypatch.setattr(eacbp.evidence.snapshot, "_verify_registry_scope", lambda *args: order.append("verified"))
    # A missing reference directory must never be reached, even with a copied
    # success summary. Native snapshot/hash/audit admission is the actual gate.
    with pytest.raises(ValueError, match="snapshot|native task|audit admission"):
        evaluate("pbmc3k", run, tmp_path / "answers-not-present", expected_image="sha256:" + "a" * 64)
    assert order == ["verified"]


@pytest.fixture
def execution_case(tmp_path):
    root = tmp_path / "curriculum"
    run = root / "runs" / "native-run"
    reference = root / "sealed_answers"
    for path in (run, reference, root / "datasets", root / "cache",
                 root / "dependencies/site-packages", root / "dependencies/linux-site-packages"):
        path.mkdir(parents=True, exist_ok=True)
    raw = root / "datasets" / "public.h5ad"
    raw.write_bytes(b"public integer counts")
    digest = sha256(raw)
    (run / "source.json").write_text(json.dumps({"path": "/input/data.h5ad", "sha256": digest,
                                               "sha256_before_import": digest, "sha256_after_import": digest}), encoding="utf-8")
    snapshot = SimpleNamespace(config={"threads": 2, "autonomous_analysis": {"question": "Explore the data", "source_sha256": digest}},
                               manifest=SimpleNamespace(biological_design=SimpleNamespace(species="mouse", tissue="pancreas"),
                                                        constraints=SimpleNamespace(max_runtime_hours=12.)))
    image = "sha256:" + "a" * 64
    command = ["docker", "run", "--rm", "--network", "none", "--read-only", "--user", "10001:10001"]
    for source, target, readonly in ((raw, "/input/data.h5ad", True), (root / "runs", "/results", False),
                                     (root / "cache", "/scratch", False), (root / "dependencies/site-packages", "/deps", True),
                                     (root / "dependencies/linux-site-packages", "/linux-deps", True)):
        command.extend(["--mount", f"type=bind,source={source.as_posix()},target={target}" + (",readonly" if readonly else "")])
    command += ["--env", "PYTHONPATH=/linux-deps:/deps", "--entrypoint", "/opt/venv/bin/eacbp", image,
                "analyze", "--question", "Explore the data", "--data", "/input/data.h5ad", "--species", "mouse",
                "--tissue", "pancreas", "--run-dir", "/results/native-run", "--threads", "2"]
    receipt_path = run.parent / f"{run.name}-execution-command.json"
    receipt = {"project": "pancreas_velocity", "command": command, "answer_reference_mounted": False}
    def save():
        receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    save()
    return SimpleNamespace(root=root, run=run, reference=reference, raw=raw, snapshot=snapshot,
                           image=image, receipt=receipt, save=save)


def test_exact_frozen_runner_receipt_is_accepted(execution_case):
    case = execution_case
    result = _execution_boundary("pancreas_velocity", case.run, case.reference, case.snapshot, case.image)
    assert result["expected_image"] == result["receipt_image"] == case.image
    assert "not a signed" in result["attestation_scope"]


@pytest.mark.parametrize("fault", ["run", "image", "project", "source", "network", "unknown", "missing_image", "runtime"])
def test_changed_runner_identity_or_unknown_shapes_are_rejected(execution_case, fault):
    case = execution_case
    command = case.receipt["command"]
    expected_image = case.image
    if fault == "run":
        command[command.index("--run-dir") + 1] = "/results/a-different-run"
    elif fault == "image":
        command[command.index(case.image)] = "sha256:" + "b" * 64
    elif fault == "project":
        case.receipt["project"] = "pbmc3k"
    elif fault == "source":
        case.raw.write_bytes(b"replaced counts")
    elif fault == "network":
        command[command.index("--network") + 1] = "host"
    elif fault == "unknown":
        command[2:2] = ["--privileged"]
    elif fault == "runtime":
        command.extend(["--max-runtime-hours", "999"])
    else:
        expected_image = None
    case.save()
    with pytest.raises(ValueError):
        _execution_boundary("pancreas_velocity", case.run, case.reference, case.snapshot, expected_image)


@pytest.mark.parametrize("style", ["mount", "volume"])
def test_resolved_reference_alias_is_rejected_for_every_mount_syntax(execution_case, monkeypatch, style):
    case = execution_case
    # Model the host resolver's symlink/junction result without requiring the
    # privilege to create a Windows symlink in a CI account.
    alias = case.root / "innocent-name"
    original_resolve = Path.resolve
    def resolve(path, *args, **kwargs):
        if path == alias:
            return case.reference
        return original_resolve(path, *args, **kwargs)
    monkeypatch.setattr(Path, "resolve", resolve)
    if style == "mount":
        case.receipt["command"][2:2] = ["--mount", f"type=bind,source={alias.as_posix()},target=/hidden,readonly"]
    else:
        case.receipt["command"][2:2] = ["--volume", f"{alias.as_posix()}:/hidden:ro"]
    case.save()
    with pytest.raises(ValueError, match="reference answers through a path or alias"):
        _execution_boundary("pancreas_velocity", case.run, case.reference, case.snapshot, case.image)


def test_reference_parent_mount_is_rejected_without_name_heuristics(execution_case):
    case = execution_case
    case.receipt["command"][2:2] = ["-v", f"{case.root.as_posix()}:/everything:ro"]
    case.save()
    with pytest.raises(ValueError, match="reference answers"):
        _execution_boundary("pancreas_velocity", case.run, case.reference, case.snapshot, case.image)


def test_supported_volume_binds_match_the_mount_form(execution_case):
    case = execution_case
    command = case.receipt["command"]
    for i, item in enumerate(command):
        if item == "--mount":
            fields = command[i + 1].split(",")
            values = dict(field.split("=", 1) for field in fields if "=" in field)
            command[i] = "--volume"
            command[i + 1] = f"{values['source']}:{values['target']}:" + ("ro" if "readonly" in fields else "rw")
    case.save()
    _execution_boundary("pancreas_velocity", case.run, case.reference, case.snapshot, case.image)


@pytest.mark.parametrize("arguments", [
    ["--volume", "opaque-named-volume:/hidden:ro"],
    ["--mount", "type=volume,source=opaque,target=/hidden"],
    ["--mount", "type=bind,source=/missing,target=/hidden,unknown=true"],
    ["--mount", "type=bind,source=/missing,target=/hidden,readonly=unknown"],
    ["--mount=type=bind,source=/missing,target=/hidden"],
])
def test_unresolvable_or_unsupported_mount_shapes_fail_closed(execution_case, arguments):
    case = execution_case
    case.receipt["command"][2:2] = arguments
    case.save()
    with pytest.raises((ValueError, FileNotFoundError)):
        _execution_boundary("pancreas_velocity", case.run, case.reference, case.snapshot, case.image)


def test_raw_count_identity_accepts_sparse_ordered_subset_and_rejects_foreign_ids():
    matrix = sparse.csr_matrix([[1, 2], [3, 4]])
    source = SimpleNamespace(X=matrix, obs_names=pd.Index(["a", "b"]), var_names=pd.Index(["g1", "g2"]))
    data = SimpleNamespace(obs=pd.DataFrame(index=["b"]), var=pd.DataFrame(index=["g2"]),
                           layers={"counts": sparse.csr_matrix([[4]])})
    assert _source_subset(source, data)
    data.layers["counts"] = sparse.csr_matrix([[40]])
    assert not _source_subset(source, data)
    data.obs.index = ["foreign"]
    assert not _source_subset(source, data)
    assert counts_equal(matrix, matrix.toarray())


def test_independent_donor_effects_keep_discordant_donors_and_missing_coverage():
    obs = pd.DataFrame({"population": ["x"] * 4, "donor": ["d1", "d1", "d2", "d2"],
                        "condition": ["stim", "ctrl", "stim", "ctrl"]}, index=["a", "b", "c", "d"])
    source = SimpleNamespace(obs=obs, obs_names=obs.index, var_names=pd.Index(["response", "background"]),
                             X=sparse.csr_matrix([[10, 100], [1, 100], [1, 100], [10, 100]]))
    spec = {"label_key": "population", "donor_key": "donor", "condition_key": "condition",
            "tested_condition": "stim", "reference_condition": "ctrl"}
    result = _independent_donor_effects(source, spec, "x", ["response"], obs.index)
    assert result["complete_donor_coverage"]
    assert result["genes"]["response"]["positive_donor_fraction"] == .5
    assert np.isclose(result["genes"]["response"]["mean_log_ratio"], 0)
    incomplete = _independent_donor_effects(source, spec, "x", ["response"], ["a", "b", "c"])
    assert not incomplete["complete_donor_coverage"]


def test_leaveout_primary_and_summary_coverage_are_separately_required():
    donors, genes = {"d1", "d2"}, {"g1", "g2"}
    folds = pd.DataFrame([{"left_out_donor": donor, "gene": gene, "removed_sample_count": 0 if donor == "__full_model__" else 2,
                           "n_pseudobulk_samples": 4 if donor == "__full_model__" else 2}
                          for donor in ["__full_model__", "d1", "d2"] for gene in sorted(genes)])
    summary = pd.DataFrame({"gene": sorted(genes), "n_requested_fits": 2, "n_successful_fits": 2,
                            "all_requested_fits_completed": True, "complete_fit_coverage": 1.})
    good = leaveout_integrity(folds, summary, donors, genes, paired=True)
    assert all(row["passed"] for row in good["checks"])
    missing_fold = leaveout_integrity(folds.loc[~folds["left_out_donor"].eq("d2")], summary, donors, genes, paired=True)
    assert not missing_fold["checks"][0]["passed"]
    incomplete_summary = summary.copy()
    incomplete_summary["n_successful_fits"] = 1
    incomplete_summary["complete_fit_coverage"] = .5
    incomplete_summary["all_requested_fits_completed"] = False
    result = leaveout_integrity(folds, incomplete_summary, donors, genes, paired=True)
    assert result["checks"][0]["passed"]
    assert not result["checks"][1]["passed"]
    # The summary cannot silently stand in for primary donor-by-gene folds.
    mistaken_primary = leaveout_integrity(summary, summary, donors, genes, paired=True)
    assert not mistaken_primary["checks"][0]["passed"]
