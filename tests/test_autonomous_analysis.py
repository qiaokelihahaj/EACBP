"""Native analyze contracts: observed design, answer isolation and lifecycle use."""

from copy import deepcopy
import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp import cli
from eacbp.application import autonomous, study_service
from eacbp.capabilities.sc_data import SCData


def data(obs=None):
    if obs is None:
        obs = pd.DataFrame(index=[f"c{i}" for i in range(8)])
    return SCData(sparse.csr_matrix(np.arange(len(obs) * 4).reshape(len(obs), 4) % 5), obs,
                  pd.DataFrame(index=["MT-A", "B", "C", "D"]))


def test_exploration_quarantines_labels_embeddings_roots_and_hvg():
    payload = data()
    payload.obs["cell_type"] = "reference_answer"
    payload.obs["paul15_clusters"] = "7MEP"
    payload.obs["louvain"] = "answer"
    payload.obs["percent_mito"] = 99.
    payload.var["highly_variable"] = True
    payload.obsm["X_pca"] = np.ones((8, 2))
    payload.obsp["connectivities"] = sparse.eye(8)
    payload.uns.update(iroot=5, paga={"answer": True}, neighbors={"answer": True})
    settings = {"question": "Discover populations and developmental trajectories", "observations": autonomous.observe_input(payload, question="Discover populations and developmental trajectories")}
    prepared, receipt = autonomous.prepare_input(payload, settings)
    assert {"cell_type", "paul15_clusters", "louvain", "percent_mito"}.isdisjoint(prepared.obs)
    assert "highly_variable" not in prepared.var
    assert not prepared.obsm and not prepared.obsp
    assert {"iroot", "paga", "neighbors"}.isdisjoint(prepared.uns)
    assert {"obs.cell_type", "obs.paul15_clusters", "var.highly_variable", "uns.iroot", "obsm.X_pca"} <= set(receipt["removed_analysis_fields"])
    assert receipt["cell_type_source"] is None
    assert "cell_type" in payload.obs and "X_pca" in payload.obsm
    np.testing.assert_array_equal(prepared.layers["counts"].toarray(), payload.X.toarray())


def test_observed_groups_are_retained_only_with_explicit_question_context():
    payload = data(pd.DataFrame({"celltype": ["T"] * 8}, index=[f"c{i}" for i in range(8)]))
    prepared, receipt = autonomous.prepare_input(payload, {"question": "Compare response across existing observed cell populations"})
    assert prepared.obs["cell_type"].tolist() == ["T"] * 8
    assert receipt["cell_type_source"] == "celltype"
    prepared, _ = autonomous.prepare_input(payload, {"question": "Analyze unannotated cell types"})
    assert "celltype" not in prepared.obs and "cell_type" not in prepared.obs


def test_generic_analysis_exports_are_quarantined_without_changing_counts_or_design():
    payload = data(pd.DataFrame({"condition": ["ctrl", "stim"] * 4,
                                 "donor": ["d1", "d1", "d2", "d2", "d3", "d3", "d4", "d4"],
                                 "batch": ["library1"] * 8}, index=[f"c{i}" for i in range(8)]))
    poisoned_obs = ("tsne1", "tsne2", "UMAP_1", "PC1", "integrated_snn_res.0.4", "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "root_cells", "end_points", "root_scores", "end_scores")
    for name in poisoned_obs:
        payload.obs[name] = 999.
    payload.var["highly_variable_genes"] = [False, True, False, True]
    payload.var["gene_ids"] = [f"id{i}" for i in range(payload.n_vars)]
    payload.uns["clusters_colors"] = ["reference_answer"]
    payload.uns["clusters_coarse_colors"] = ["reference_answer"]
    payload.uns["accession"] = "objective_source_identifier"
    original = payload.X.toarray().copy()
    prepared, receipt = autonomous.prepare_input(payload, {"question": "Discover cell populations"})
    assert set(poisoned_obs).isdisjoint(prepared.obs)
    assert "highly_variable_genes" not in prepared.var
    assert {"clusters_colors", "clusters_coarse_colors"}.isdisjoint(prepared.uns)
    assert {f"obs.{name}" for name in poisoned_obs} <= set(receipt["removed_analysis_fields"])
    assert {"var.highly_variable_genes", "uns.clusters_colors", "uns.clusters_coarse_colors"} <= set(receipt["removed_analysis_fields"])
    np.testing.assert_array_equal(prepared.X.toarray(), original)
    np.testing.assert_array_equal(prepared.layers["counts"].toarray(), original)
    for name in ("condition", "donor", "batch"):
        assert prepared.obs[name].tolist() == payload.obs[name].tolist()
    assert prepared.var["gene_ids"].tolist() == payload.var["gene_ids"].tolist()
    assert prepared.uns["accession"] == "objective_source_identifier"


def test_alias_design_detects_all_eligible_targets_and_real_pairing():
    obs = pd.DataFrame([
        {"cell_type": target, "label": condition, "mouse_id": donor}
        for target, donors in (("T", ("d1", "d2")), ("B", ("d1",)), ("NK", ("d1", "d2", "d3")))
        for condition in ("stim", "ctrl") for donor in donors
    ])
    observation = autonomous.observe_input(data(obs), question="Compare IFN response in existing observed cell populations")
    assert observation["obs_column_map"]["condition"] == "label"
    assert observation["obs_column_map"]["donor"] == "mouse_id"
    assert observation["contrast"]["tested"] == "stim"
    targets = {item["cell_type"]: item for item in observation["targets"]}
    assert targets["T"]["eligible"] and targets["NK"]["eligible"]
    assert targets["T"]["paired"] and not targets["B"]["eligible"]
    reordered = obs.iloc[::-1].copy()
    assert autonomous.observe_input(data(reordered), question="Compare IFN response in existing observed cell populations")["contrast"] == observation["contrast"]


def test_unmatched_paired_donors_never_become_independent_samples():
    obs = pd.DataFrame({"cell_type": ["T"] * 4, "condition": ["ctrl", "ctrl", "stim", "stim"], "donor": ["d1", "d2", "d1", "d3"]})
    target = autonomous.observe_input(data(obs), question="Compare observed cell types")["targets"][0]
    assert target["paired"] and not target["eligible"]


def test_named_counts_are_authoritative_and_transformed_x_is_rejected():
    payload = data()
    payload.X = payload.X.astype(float)
    payload.X.data += .2
    with pytest.raises(ValueError, match="Unsupported expression representation"):
        autonomous.observe_input(payload, question="Discover cell populations")
    payload.layers["counts"] = sparse.csr_matrix(np.ones(payload.shape, dtype=int))
    assert autonomous.observe_input(payload, question="Discover cell populations")["counts_source"] == "layers.counts"
    payload.layers["counts"].data = payload.layers["counts"].data.astype(float) + .5
    with pytest.raises(ValueError, match="Named counts layer is invalid"):
        autonomous.observe_input(payload, question="Discover cell populations")


def test_metadata_ambiguity_and_changed_observations_fail_closed():
    payload = data(pd.DataFrame({"condition": ["ctrl"] * 4 + ["stim"] * 4, "treatment": ["x"] * 8}))
    with pytest.raises(ValueError, match="Ambiguous condition"):
        autonomous.observe_input(payload, question="Compare observed cell types")
    payload = data()
    observation = autonomous.observe_input(payload, question="Discover populations")
    changed = deepcopy(observation)
    changed["n_cells"] += 1
    with pytest.raises(ValueError, match="observations changed"):
        autonomous.prepare_input(payload, {"question": "Discover populations", "observations": changed})


def test_native_analysis_calls_shared_lifecycle_and_never_accepts_config(tmp_path, monkeypatch):
    source = tmp_path / "input.h5ad"
    source.write_bytes(b"test file identity")
    obs = pd.DataFrame({"cell_type": ["T"] * 8, "label": ["stim", "ctrl"] * 4, "mouse_id": ["d1", "d1", "d2", "d2", "d3", "d3", "d4", "d4"]})
    payload = data(obs)
    monkeypatch.setattr(study_service, "_import_h5ad", lambda path: payload)
    calls = []
    monkeypatch.setattr(study_service, "run_study", lambda **kwargs: calls.append(kwargs) or {"status": "success"})
    result = autonomous.analyze_study(question="Compare response in existing observed cell populations", data=source, species="human", tissue="PBMC", run_dir=tmp_path / "run")
    assert result["status"] == "success" and len(calls) == 1
    assert result == {"status": "success"}
    call = calls[0]
    assert call["data"] == source and call["manifest"].biological_design.target_cell_types == ["T"]
    settings = call["config"]["autonomous_analysis"]
    assert settings["intent"] == "condition_comparison"
    assert settings["source_sha256"] == study_service._hash_file(source)
    assert call["config"]["capability_parameters"]["deg"]["paired"] is True
    assert call["config"]["analysis_extensions"]["donor_sensitivity"] is True
    from eacbp.capabilities import create_default_capability_registry
    from eacbp.orchestrator.planning import build_study_tasks
    tasks = build_study_tasks(call["manifest"], call["config"], create_default_capability_registry())
    scoped = [task for task in tasks if task.capability in {"subset_cells", "deg", "donor_sensitivity"}]
    assert len(scoped) == 3 and all(task.parameters["target_cell_type"] == "T" for task in scoped)
    with pytest.raises(TypeError, match="config"):
        autonomous.analyze_study(question="Discover populations", data=source, species="human", tissue="PBMC", run_dir=tmp_path / "other", config={"marker_reference": {"answer": ["B"]}})


def test_analyze_cli_excludes_answer_and_method_inputs():
    args = ["analyze", "--question", "Discover populations", "--data", "raw.h5ad", "--species", "human", "--tissue", "PBMC", "--run-dir", "run"]
    assert cli.build_parser().parse_args(args).command == "analyze"
    for option in ("--config", "--manifest", "--root-cell", "--targets", "--reference", "--min-genes", "--method"):
        with pytest.raises(SystemExit):
            cli.build_parser().parse_args([*args, option, "injected-answer"])


def test_condition_question_without_group_context_records_gap_and_preserves_estimand(tmp_path, monkeypatch):
    source = tmp_path / "input.h5ad"
    source.write_bytes(b"test file identity")
    obs = pd.DataFrame({"cell_type": ["T"] * 8, "condition": ["case", "ctrl"] * 4,
                        "donor": ["d1", "d1", "d2", "d2", "d3", "d3", "d4", "d4"]})
    monkeypatch.setattr(study_service, "_import_h5ad", lambda path: data(obs))
    calls = []
    monkeypatch.setattr(study_service, "run_study", lambda **kwargs: calls.append(kwargs) or {"status": "success"})
    result = autonomous.analyze_study(question="Compare treatment response across cells", data=source, species="human", tissue="blood", run_dir=tmp_path / "run")
    config = calls[0]["config"]
    assert config["autonomous_analysis"]["intent"] == "exploration"
    assert "mixed-cell replacement estimand" in config["autonomous_analysis"]["scientific_gaps"][0]
    assert calls[0]["manifest"].biological_design.target_cell_types == []
    assert "advanced_analysis" not in config
    assert result == {"status": "success"}


def test_batch_condition_confounding_disables_correction(tmp_path, monkeypatch):
    source = tmp_path / "input.h5ad"
    source.write_bytes(b"test file identity")
    obs = pd.DataFrame({"condition": ["stim"] * 4 + ["ctrl"] * 4,
                        "batch": ["b1"] * 4 + ["b2"] * 4})
    monkeypatch.setattr(study_service, "_import_h5ad", lambda path: data(obs))
    calls = []
    monkeypatch.setattr(study_service, "run_study", lambda **kwargs: calls.append(kwargs) or {"status": "success"})
    autonomous.analyze_study(question="Discover populations", data=source, species="human", tissue="blood", run_dir=tmp_path / "run")
    config = calls[0]["config"]
    assert config["autonomous_analysis"]["observations"]["batch_condition_confounded"] is True
    assert config["method_overrides"]["integration"] == "no_correction_v1"


def test_no_condition_never_generates_a_differential_contrast(tmp_path, monkeypatch):
    source = tmp_path / "input.h5ad"
    source.write_bytes(b"test file identity")
    monkeypatch.setattr(study_service, "_import_h5ad", lambda path: data())
    calls = []
    monkeypatch.setattr(study_service, "run_study", lambda **kwargs: calls.append(kwargs) or {"status": "success"})
    autonomous.analyze_study(question="Discover cell populations and developmental trajectories", data=source, species="mouse", tissue="marrow", run_dir=tmp_path / "run")
    config = calls[0]["config"]
    receipt = config["autonomous_analysis"]
    assert receipt["observations"]["contrast"] is None and receipt["intent"] == "exploration"
    assert "deg" not in config["capability_parameters"] and "advanced_analysis" not in config
    assert config["analysis_extensions"] == {"cluster_markers": True, "qc_review": True}
    assert any("DPT" in gap for gap in receipt["scientific_gaps"])


def test_velocity_uses_validated_kinetic_layers_and_quarantines_prior_outputs(tmp_path, monkeypatch):
    source = tmp_path / "input.h5ad"
    source.write_bytes(b"test file identity")
    payload = data()
    payload.layers["spliced"] = payload.X.copy()
    payload.layers["unspliced"] = payload.X.copy()
    payload.layers["velocity"] = payload.X.astype(float)
    payload.obs["clusters"] = "reference"
    payload.obs["S_score"] = .9
    payload.obs["G2M_score"] = .8
    payload.obs["latent_time"] = .7
    payload.uns["iroot"] = 2
    payload.uns["velocity_params"] = {"method": "precomputed"}
    question = "Recover endocrine populations and infer RNA velocity developmental direction"
    monkeypatch.setattr(study_service, "_import_h5ad", lambda path: payload)
    calls = []
    monkeypatch.setattr(study_service, "run_study", lambda **kwargs: calls.append(kwargs) or {"status": "success"})
    result = autonomous.analyze_study(question=question, data=source, species="mouse", tissue="pancreas", run_dir=tmp_path / "run")
    config = calls[0]["config"]
    assert config["analysis_extensions"]["rna_velocity"] is True
    assert result == {"status": "success"}
    assert config["autonomous_analysis"]["limitations"] and not config["autonomous_analysis"]["scientific_gaps"]
    assert config["autonomous_analysis"]["observations"]["rna_velocity"]["available"] is True
    prepared, receipt = autonomous.prepare_input(payload, config["autonomous_analysis"])
    assert set(prepared.layers) == {"counts", "spliced", "unspliced"}
    assert {"clusters", "S_score", "G2M_score", "latent_time"}.isdisjoint(prepared.obs)
    assert "iroot" not in prepared.uns and "velocity_params" not in prepared.uns
    assert receipt["retained_kinetic_layers"] == ["spliced", "unspliced"]
    assert "layers.velocity" in receipt["removed_analysis_fields"]
    np.testing.assert_array_equal(prepared.layers["unspliced"].toarray(), payload.layers["unspliced"].toarray())


def test_missing_or_invalid_kinetic_layers_make_the_velocity_request_partial(tmp_path, monkeypatch):
    source = tmp_path / "input.h5ad"
    source.write_bytes(b"test file identity")
    payload = data()
    payload.layers["spliced"] = payload.X.copy()
    payload.layers["unspliced"] = payload.X.astype(float)
    payload.layers["unspliced"].data += .5
    monkeypatch.setattr(study_service, "_import_h5ad", lambda path: payload)
    calls = []
    monkeypatch.setattr(study_service, "run_study", lambda **kwargs: calls.append(kwargs) or {"status": "success"})
    result = autonomous.analyze_study(question="Infer RNA velocity developmental direction", data=source, species="mouse", tissue="pancreas", run_dir=tmp_path / "run")
    assert "rna_velocity" not in calls[0]["config"]["analysis_extensions"]
    assert result == {"status": "success"}
    assert any("RNA velocity is not estimable" in gap for gap in calls[0]["config"]["autonomous_analysis"]["scientific_gaps"])
    prepared, receipt = autonomous.prepare_input(payload, calls[0]["config"]["autonomous_analysis"])
    assert set(prepared.layers) == {"counts"} and receipt["retained_kinetic_layers"] == []


def test_requested_comparison_without_condition_design_is_partial_and_cli_fails(tmp_path, monkeypatch, capsys):
    source = tmp_path / "input.h5ad"
    source.write_bytes(b"test file identity")
    monkeypatch.setattr(study_service, "_import_h5ad", lambda path: data())
    calls = []
    monkeypatch.setattr(study_service, "run_study", lambda **kwargs: calls.append(kwargs) or {"status": "success", "analysis_request_status": "partial"})
    result = autonomous.analyze_study(question="Compare treatment response in existing observed cell populations", data=source, species="human", tissue="blood", run_dir=tmp_path / "run")
    assert result["analysis_request_status"] == "partial"
    assert any("two-condition experimental design" in gap for gap in calls[0]["config"]["autonomous_analysis"]["scientific_gaps"])
    monkeypatch.setattr(autonomous, "analyze_study", lambda **kwargs: result)
    assert cli.main(["analyze", "--question", "Compare treatment response", "--data", str(source), "--species", "human", "--tissue", "blood", "--run-dir", str(tmp_path / "cli")]) == 1
    assert '"analysis_request_status": "partial"' in capsys.readouterr().out


def test_request_coverage_requires_all_population_fits_and_durable_admission():
    from eacbp.auditor.base import ValidationReport
    from eacbp.schemas.audit import stable_fingerprint
    from eacbp.schemas.study import BiologicalDesign, StudyManifest
    from eacbp.schemas.task import ScientificResult, ScientificResultStatus, TaskResult, TaskStatus

    manifest = StudyManifest(study_id="coverage", biological_design=BiologicalDesign(species="human", tissue="blood", target_cell_types=["T", "B"]))
    config = {"autonomous_analysis": {"intent": "condition_comparison", "question": "Compare observed cell populations", "scientific_gaps": []},
              "analysis_extensions": {"cluster_markers": True}}
    history, reports, records = [], [], []

    def admitted(capability, target=None):
        task_id = capability + ("_" + target if target else "")
        params = {"cell_type": target} if capability == "subset_cells" else {"target_cell_type": target} if target else {}
        output = f"table://coverage/{task_id}/v1"
        inputs = [f"table://coverage/subset_cells_{target}/v1"] if capability in {"deg", "donor_sensitivity"} else []
        result = TaskResult(task_id=task_id, capability=capability, method_used="real_method", status=TaskStatus.SUCCESS,
                            input_artifacts=inputs, output_artifacts=[output],
                            metrics={"all_requested_fits_completed": True} if capability == "donor_sensitivity" else {})
        if capability in {"deg", "donor_sensitivity"}:
            result.scientific_result = ScientificResult(status=ScientificResultStatus.ESTIMATED_INCONCLUSIVE,
                                                       summary="Estimates produced with no supported features", n_features_tested=3, n_features_supported=0)
        history.append(result)
        reports.append(ValidationReport(auditor_name="independent", target_task_id=task_id, overall_passed=True))
        context = {"capability": capability, "method": "real_method", "parameters": params}
        records.append(SimpleNamespace(task_id=task_id, overall_passed=True, stop_rule_triggered=False,
                                       output_artifacts=[output], contract=context, contract_fingerprint=stable_fingerprint(context)))

    for capability in ("dataset_audit", "qc", "normalization", "integration", "clustering", "cluster_markers"):
        admitted(capability)
    for target in ("T", "B"):
        for capability in ("subset_cells", "deg", "donor_sensitivity"):
            admitted(capability, target)
    registry = SimpleNamespace(list_audit_records=lambda **kwargs: records)
    orchestrator = SimpleNamespace(artifact_registry=registry, task_history=history, audit_reports=reports)
    certify = lambda: study_service._autonomous_request_summary({"status": "success"}, manifest=manifest, config=config, orchestrator=orchestrator)
    assert certify()["analysis_request_status"] == "complete"
    from eacbp.capabilities import create_default_capability_registry
    orchestrator.capability_registry = create_default_capability_registry()
    normalization = next(result for result in history if result.capability == "normalization")
    normalization_record = next(record for record in records if record.contract["capability"] == "normalization")
    normalization.method_used = "library_size_log1p_v1"
    normalization_record.contract["method"] = "sc_normalize_log1p_v1"
    normalization_record.contract_fingerprint = stable_fingerprint(normalization_record.contract)
    assert certify()["analysis_request_status"] == "complete"
    normalization_record.contract["method"] = "unregistered_normalization_alias"
    assert certify()["analysis_request_status"] == "partial"
    normalization_record.contract["method"] = "sc_normalize_log1p_v1"
    integration = next(result for result in history if result.capability == "integration")
    integration_record = next(record for record in records if record.contract["capability"] == "integration")
    integration.method_used = "no_correction_v1"
    integration_record.contract["method"] = "no_correction_v1"
    integration_record.contract_fingerprint = stable_fingerprint(integration_record.contract)
    assert certify()["analysis_request_status"] == "complete"
    # Even another registered alias for the same method cannot rewrite the
    # saved receipt context after independent admission.
    integration_record.contract["method"] = "no_correction"
    assert certify()["analysis_request_status"] == "partial"
    integration_record.contract["method"] = "harmonypy_v1"
    integration_record.contract_fingerprint = stable_fingerprint(integration_record.contract)
    assert certify()["analysis_request_status"] == "partial"
    integration_record.contract["method"] = "no_correction_v1"
    integration_record.contract_fingerprint = stable_fingerprint(integration_record.contract)
    history[-1].metrics["all_requested_fits_completed"] = False
    assert certify()["analysis_request_status"] == "partial"
    history[-1].metrics["all_requested_fits_completed"] = True
    records.pop()
    failed = certify()
    assert failed["analysis_request_status"] == "partial"
    assert any("donor_sensitivity for population 'B'" in gap for gap in failed["scientific_gaps"])


def test_partial_native_request_status_survives_snapshot_resume_and_portability(tmp_path, monkeypatch):
    import anndata as ad
    from eacbp.orchestrator import loop

    original = loop.build_study_tasks
    monkeypatch.setattr(loop, "build_study_tasks", lambda *args, **kwargs: [task for task in original(*args, **kwargs) if task.capability == "dataset_audit"])
    source = tmp_path / "source.h5ad"
    payload = data()
    ad.AnnData(payload.X, obs=payload.obs, var=payload.var).write_h5ad(source)
    run_dir = tmp_path / "run"
    result = autonomous.analyze_study(question="Discover cell populations", data=source, species="human", tissue="blood", run_dir=run_dir)
    assert result["status"] == "success" and result["analysis_request_status"] == "partial"
    assert any("cluster_markers" in gap for gap in result["scientific_gaps"])
    saved = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    snapshot = json.loads((run_dir / "snapshot.json").read_text(encoding="utf-8"))
    assert saved["analysis_completion_checks"] == snapshot["summary"]["analysis_completion_checks"]
    assert saved["requirement_completion_checks"] == snapshot["summary"]["requirement_completion_checks"]
    assert snapshot["summary"]["analysis_request_status"] == "partial"
    assert "Requested analysis: `partial`" in (run_dir / "report.md").read_text(encoding="utf-8")
    resumed = study_service.resume_study(run_dir=run_dir)
    assert resumed["status"] == "success" and resumed["analysis_request_status"] == "partial"
    assert resumed["scientific_gaps"] == saved["scientific_gaps"]
    assert resumed["requirement_completion_checks"] == saved["requirement_completion_checks"]
    bundle = tmp_path / "study.zip"
    study_service.export_study(run_dir=run_dir, output=bundle)
    imported = tmp_path / "imported"
    study_service.import_study(bundle=bundle, run_dir=imported)
    portable = json.loads((imported / "snapshot.json").read_text(encoding="utf-8"))
    portable_summary = json.loads((imported / "summary.json").read_text(encoding="utf-8"))
    assert portable["summary"]["analysis_request_status"] == portable_summary["analysis_request_status"] == "partial"
    assert portable_summary["requirement_completion_checks"] == saved["requirement_completion_checks"]
    assert json.loads((imported / "autonomy.json").read_text(encoding="utf-8"))["question"] == "Discover cell populations"
    assert (imported / "input_preparation.json").is_file()
