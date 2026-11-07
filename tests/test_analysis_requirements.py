"""Complex question scope must survive eligibility and partial execution."""

from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.application import autonomous, study_service
from eacbp.application.analysis_requirements import (
    evaluate_requirement_coverage, freeze_question_requirements, requested_dimensions,
)
from eacbp.capabilities.sc_data import SCData


def _data(obs):
    return SCData(sparse.csr_matrix(np.arange(len(obs) * 4).reshape(len(obs), 4) % 5),
                  obs, pd.DataFrame(index=["MT-A", "B", "C", "D"]))


def _observation(*, conditions=("ctrl", "stim"), underpowered=False):
    obs = pd.DataFrame([
        {"cell_type": target, "condition": condition, "donor": donor}
        for target in ("B", "T") for condition in conditions
        for donor in (("d1",) if underpowered and target == "B" else ("d1", "d2", "d3"))
    ])
    return autonomous.observe_input(_data(obs), question="Compare existing observed cell populations")


def _admitted(capability, *, target=None, contrast=None, scientific_status="estimated_inconclusive", network_kind=None):
    task_id = capability + ("_" + target if target else "") + ("_" + str(contrast["tested"]) if contrast else "")
    parameters = {"cell_type": target} if capability == "subset_cells" else {"target_cell_type": target} if target else {}
    if contrast:
        parameters.update(condition_a=contrast["tested"], condition_b=contrast["reference"])
    if network_kind:
        parameters["network_kind"] = network_kind
    metrics = {"all_requested_fits_completed": True} if capability == "donor_sensitivity" else {}
    scientific = SimpleNamespace(status=scientific_status, n_features_tested=3)
    return (SimpleNamespace(task_id=task_id, capability=capability, status="success", method_used="real_method",
                            output_artifacts=[f"table://study/{task_id}/v1"],
                            input_artifacts=[f"table://study/subset_cells_{target}/v1"] if target and capability != "subset_cells" else [],
                            metrics=metrics, scientific_result=scientific),
            {"capability": capability, "method": "real_method", "parameters": parameters})


def _native_results(contract, *, omit_target=None):
    records = [_admitted(capability) for capability in
               ("dataset_audit", "qc", "normalization", "integration", "clustering", "cluster_markers")]
    for target in contract["population_universe"]:
        if target == omit_target:
            continue
        records.append(_admitted("subset_cells", target=target))
        for contrast in contract["contrast_universe"]:
            if contrast["resolved"]:
                records.extend(_admitted(capability, target=target, contrast=contrast)
                               for capability in ("deg", "donor_sensitivity"))
    return records


def test_frozen_universe_retains_ineligible_population_and_never_claims_complete():
    observation = _observation(underpowered=True)
    contract = freeze_question_requirements(question="Compare existing observed cell populations", observations=observation)
    assert contract["population_universe"] == ["B", "T"]
    assert {requirement["target_cell_type"] for requirement in contract["requirements"]
            if requirement["dimension"] == "population_expression"} == {"B", "T"}
    coverage = evaluate_requirement_coverage(contract, _native_results(contract, omit_target="B"))
    assert coverage["analysis_request_status"] == "partial"
    assert any(check["target_cell_type"] == "B" and check["status"] == "not_estimable"
               for check in coverage["requirement_completion_checks"])
    assert any("'B'" in gap for gap in coverage["scientific_gaps"])


def test_missing_one_population_remains_partial_despite_successful_other_population():
    contract = freeze_question_requirements(question="Compare existing observed cell populations", observations=_observation())
    coverage = evaluate_requirement_coverage(contract, _native_results(contract, omit_target="B"))
    assert coverage["analysis_request_status"] == "partial"
    assert any(check["target_cell_type"] == "B" and check["status"] == "missing_or_unadmitted"
               for check in coverage["requirement_completion_checks"])


def test_audited_null_estimates_can_complete_original_binary_request():
    contract = freeze_question_requirements(question="Compare existing observed cell populations", observations=_observation())
    coverage = evaluate_requirement_coverage(contract, _native_results(contract))
    assert coverage["analysis_request_status"] == "complete"
    assert "estimated_inconclusive" in {status for check in coverage["requirement_completion_checks"]
                                       for status in check["scientific_statuses"]}


@pytest.mark.parametrize("status", ["not_estimable", "assumptions_failed"])
def test_scientific_failures_cannot_borrow_task_success(status):
    contract = freeze_question_requirements(question="Compare existing observed cell populations", observations=_observation())
    records = _native_results(contract)
    next(result for result, _ in records if result.capability == "deg").scientific_result.status = status
    coverage = evaluate_requirement_coverage(contract, records)
    assert coverage["analysis_request_status"] == "partial"
    assert status in {check["status"] for check in coverage["requirement_completion_checks"]}


def test_correct_declared_population_without_subset_lineage_is_not_completion():
    contract = freeze_question_requirements(question="Compare existing observed cell populations", observations=_observation())
    records = _native_results(contract)
    next(result for result, _ in records if result.capability == "deg").input_artifacts = ["adata://study/all_cells/v1"]
    assert evaluate_requirement_coverage(contract, records)["analysis_request_status"] == "partial"


def test_wrong_contrast_and_partial_leaveout_are_unmet():
    contract = freeze_question_requirements(question="Compare existing observed cell populations", observations=_observation())
    records = _native_results(contract)
    next(recorded for result, recorded in records if result.capability == "deg")["parameters"]["condition_a"] = "unobserved"
    next(result for result, _ in records if result.capability == "donor_sensitivity").metrics["all_requested_fits_completed"] = False
    assert evaluate_requirement_coverage(contract, records)["analysis_request_status"] == "partial"


def test_compound_complex_request_freezes_every_dimension_before_method_selection():
    observation = _observation(conditions=("ctrl", "drug1", "drug2"))
    observation.update(spatial={"available": True}, sections=["s1", "s2"], rna_velocity={"available": True})
    question = ("Compare multiple conditions in existing observed cell populations: abundance, pathways and gene programs, "
                "TF regulators, ligand receptor communication, replicated spatial niches, branching trajectories, "
                "fate probabilities, dynamical kinetics and donor robustness")
    dimensions = requested_dimensions(question)
    assert {"comparison", "abundance", "pathway", "regulator", "communication", "spatial", "replication",
            "multi_condition", "trajectory", "branching", "fate", "dynamic_kinetics", "robustness"} <= set(dimensions)
    contract = freeze_question_requirements(question=question, observations=observation)
    assert {"abundance", "pathway", "regulator", "communication", "replicated_spatial", "multi_condition",
            "branching", "fate", "dynamic_kinetics"} <= {item["dimension"] for item in contract["requirements"]}
    coverage = evaluate_requirement_coverage(contract, _native_results(contract), available_capabilities={"deg", "donor_sensitivity", "functional_activity", "liana_communication"})
    assert coverage["analysis_request_status"] == "partial"
    assert any(check["status"] == "unsupported" for check in coverage["requirement_completion_checks"])


def test_ligand_receptor_interaction_is_not_a_factorial_condition_request():
    assert "communication" in requested_dimensions("Infer ligand receptor interaction")
    assert "multi_condition" not in requested_dimensions("Infer ligand receptor interaction")


def test_population_result_cannot_replace_required_shared_exploration():
    contract = freeze_question_requirements(question="Discover populations", observations=_observation())
    records = _native_results(contract)
    next(recorded for result, recorded in records if result.capability == "clustering")["parameters"]["target_cell_type"] = "T"
    assert evaluate_requirement_coverage(contract, records)["analysis_request_status"] == "partial"


def test_three_conditions_do_not_invent_contrast_or_hide_population_scopes():
    contract = freeze_question_requirements(question="Compare multiple conditions in existing observed cell populations",
                                             observations=_observation(conditions=("ctrl", "drug1", "drug2")))
    assert contract["condition_universe"] == ["ctrl", "drug1", "drug2"]
    assert contract["contrast_universe"] == [{"tested": None, "reference": None, "resolved": False,
                                             "reason": "The question does not resolve a contrast from the observed condition universe"}]
    assert {item["target_cell_type"] for item in contract["requirements"]
            if item["dimension"] == "population_expression"} == {"B", "T"}


def test_explicit_observed_all_pairwise_and_control_scopes_are_frozen():
    observation = _observation(conditions=("ctrl", "drug1", "drug2"))
    pairs = freeze_question_requirements(question="Compare all pairwise comparisons in existing observed cell populations", observations=observation)
    assert len(pairs["contrast_universe"]) == 3
    assert all(item["tested"] in observation["conditions"] and item["reference"] in observation["conditions"]
               for item in pairs["contrast_universe"])
    control = freeze_question_requirements(question="Compare each condition versus control in existing observed cell populations", observations=observation)
    assert {(item["tested"], item["reference"]) for item in control["contrast_universe"]} == {("drug1", "ctrl"), ("drug2", "ctrl")}
    assert "multi_condition" not in {item["dimension"] for item in pairs["requirements"]}
    assert "multi_condition" not in {item["dimension"] for item in control["requirements"]}


def test_factorial_question_with_two_condition_levels_does_not_require_third_level():
    contract = freeze_question_requirements(question="Compare factorial condition interaction effects in existing observed cell populations", observations=_observation())
    joint = next(item for item in contract["requirements"] if item["dimension"] == "multi_condition")
    assert joint["prerequisite_gaps"] == []
    assert evaluate_requirement_coverage(contract, _native_results(contract))["analysis_request_status"] == "partial"


def test_one_unspecified_prior_cannot_satisfy_pathway_and_regulator_requirements():
    contract = freeze_question_requirements(question="Compare pathways and TF regulators in existing observed cell populations", observations=_observation())
    records = _native_results(contract)
    records.extend(_admitted("functional_activity", target=target, contrast=contract["contrast_universe"][0], network_kind="pathway")
                   for target in contract["population_universe"])
    coverage = evaluate_requirement_coverage(contract, records)
    assert coverage["analysis_request_status"] == "partial"
    assert all(check["status"] == "passed" for check in coverage["requirement_completion_checks"] if check["dimension"] == "pathway")
    assert all(check["status"] != "passed" for check in coverage["requirement_completion_checks"] if check["dimension"] == "regulator")


def test_contract_is_deeply_frozen_and_tampering_is_rejected():
    observation = _observation()
    contract = freeze_question_requirements(question="Compare existing observed cell populations", observations=observation, source_sha256="a" * 64)
    saved = deepcopy(contract)
    observation["targets"].clear()
    observation["obs_column_map"]["donor"] = "invented"
    assert contract == saved
    contract["population_universe"].remove("B")
    with pytest.raises(ValueError, match="identity or version"):
        evaluate_requirement_coverage(contract, [])
    with pytest.raises(ValueError, match="identity or version"):
        evaluate_requirement_coverage({"schema_version": 0, "requirements": []}, [])


def test_sample_and_library_ids_never_become_biological_donors():
    payload = _data(pd.DataFrame({"sample_id": ["s1", "s2"] * 4, "library_id": ["l1"] * 8,
                                  "condition": ["ctrl", "stim"] * 4, "cell_type": ["T"] * 8}))
    observation = autonomous.observe_input(payload, question="Compare existing observed cell populations")
    assert observation["obs_column_map"]["donor"] is None
    assert observation["obs_column_map"]["sample"] == "sample_id"
    assert observation["obs_column_map"]["library"] == "library_id"
    prepared, _ = autonomous.prepare_input(payload, {"question": "Compare existing observed cell populations"})
    assert "donor" not in prepared.obs and prepared.obs["sample"].tolist() == payload.obs["sample_id"].tolist()


def test_measured_coordinates_sections_and_time_survive_derived_answer_quarantine():
    payload = _data(pd.DataFrame({"section_id": ["s1"] * 4 + ["s2"] * 4, "timepoint": ["0", "1"] * 4,
                                  "spatial_domain": ["answer"] * 8, "velocity_pseudotime": [.3] * 8}))
    coordinates = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]] * 2)
    payload.obsm.update(spatial=coordinates, X_umap=np.ones((8, 2)))
    payload.obsp["spatial_connectivities"] = sparse.eye(8)
    payload.uns["spatial_connectivities"] = sparse.eye(8)
    observation = autonomous.observe_input(payload, question="Explore replicated spatial niches")
    prepared, receipt = autonomous.prepare_input(payload, {"question": "Explore replicated spatial niches", "observations": observation})
    np.testing.assert_array_equal(prepared.obsm["spatial"], coordinates)
    assert set(prepared.obsm) == {"spatial"} and not prepared.obsp
    assert {"spatial_domain", "velocity_pseudotime"}.isdisjoint(prepared.obs)
    assert "spatial_connectivities" not in prepared.uns
    assert prepared.obs["section"].tolist() == payload.obs["section_id"].tolist()
    assert prepared.obs["measured_time"].tolist() == payload.obs["timepoint"].tolist()
    assert receipt["retained_measured_coordinate_fields"] == ["spatial"]
    payload.obsm["spatial"][0, 0] = 50.
    with pytest.raises(ValueError, match="observations changed"):
        autonomous.prepare_input(payload, {"question": "Explore replicated spatial niches", "observations": observation})


def test_invalid_or_ambiguous_spatial_coordinates_do_not_become_available():
    payload = _data(pd.DataFrame(index=range(8)))
    payload.obsm["spatial"] = np.full((8, 2), np.nan)
    assert not autonomous.observe_input(payload, question="Explore spatial niches")["spatial"]["available"]
    payload.obsm["spatial"] = np.arange(16).reshape(8, 2)
    payload.obsm["spatial_coords"] = np.ones((8, 2))
    observation = autonomous.observe_input(payload, question="Explore spatial niches")
    assert not observation["spatial"]["available"] and "disagree" in observation["spatial"]["reason"]


def test_native_entry_freezes_scope_before_removing_underpowered_population(tmp_path, monkeypatch):
    source = tmp_path / "input.h5ad"
    source.write_bytes(b"raw input identity")
    obs = pd.DataFrame([{ "cell_type": target, "condition": condition, "donor": donor }
                        for target, donors in (("B", ("d1",)), ("T", ("d1", "d2", "d3")))
                        for condition in ("ctrl", "stim") for donor in donors])
    monkeypatch.setattr(study_service, "_import_h5ad", lambda _: _data(obs))
    calls = []
    monkeypatch.setattr(study_service, "run_study", lambda **kwargs: calls.append(kwargs) or {"status": "success"})
    autonomous.analyze_study(question="Compare pathways in existing observed cell populations", data=source,
                             species="human", tissue="blood", run_dir=tmp_path / "run")
    settings = calls[0]["config"]["autonomous_analysis"]
    assert calls[0]["manifest"].biological_design.target_cell_types == ["T"]
    assert settings["requirement_contract"]["population_universe"] == ["B", "T"]
    assert settings["requirement_contract"]["source_sha256"] == settings["source_sha256"]
    assert any("'B'" in gap for gap in settings["scientific_gaps"])
