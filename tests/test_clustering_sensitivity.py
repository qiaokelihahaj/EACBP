"""Counterexamples for frozen count-only clustering sensitivity contracts."""

from copy import deepcopy
import hashlib
import math

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.application.analysis_requirements import freeze_question_requirements
from eacbp.application.autonomous import observe_input, prepare_input
from eacbp.application.clustering_sensitivity import (
    clustering_branch_config, clustering_config_sha256,
    clustering_scientific_config_sha256, compare_clustering_sensitivity_runs,
    feature_axis_fingerprint, freeze_clustering_sensitivity_plan,
    raw_count_fingerprint, validate_clustering_branch_config,
    validate_clustering_sensitivity_plan,
)
from eacbp.capabilities.sc_data import SCData


def _fixture(*, donors=4, n_cells=80, source="layers.counts"):
    counts = sparse.csr_matrix(np.random.default_rng(19).poisson(6, (n_cells, 30)).astype(np.int64))
    counts[:min(3, n_cells), :] = 0
    counts.eliminate_zeros()
    obs = pd.DataFrame(index=[f"cell_{number}" for number in range(n_cells)])
    if donors:
        obs["donor_id"] = [f"GSM{2230757 + number % donors}" for number in range(n_cells)]
    var = pd.DataFrame(index=[f"gene_{number}" for number in range(30)])
    payload = SCData(counts, obs, var, layers={"counts": counts} if source == "layers.counts" else {})
    question = "Explore cell populations and marker genes"
    observation = observe_input(payload, question=question)
    settings = {"intent": "exploration", "question": question, "source_sha256": "a" * 64,
                "observations": observation, "policy_version": "observed_single_cell_v3",
                "requirement_contract": freeze_question_requirements(question=question, observations=observation,
                                                                       source_sha256="a" * 64),
                "input_preparation": {key: observation[key] for key in ("counts_source", "gene_name_source", "obs_column_map")}}
    config = {"mode": "real", "method_profile": "standard", "threads": 1,
              "autonomous_analysis": settings, "analysis_extensions": {"cluster_markers": True, "qc_review": True},
              "method_overrides": {"integration": "no_correction_v1"},
              "capability_parameters": {"qc": {key: observation["qc"][key] for key in ("min_genes", "max_mito_pct")},
                                        "normalization": {"n_top_genes": 20},
                                        "integration": {"n_components": 10, "random_seed": 42},
                                        "clustering": {"resolution": 1., "random_seed": 42, "n_neighbors": 15}}}
    plan = freeze_clustering_sensitivity_plan(payload=payload, observations=observation,
                                             source_sha256="a" * 64, baseline_config=config)
    return payload, observation, config, plan


def _spec(plan, branch_id):
    return next(item for item in [plan["primary"], *plan["branches"]] if item["branch_id"] == branch_id)


def _receipt(plan, config, branch_id):
    spec = _spec(plan, branch_id)
    source = plan["source_sha256"] if spec["source_kind"] == "original" else hashlib.sha256(branch_id.encode()).hexdigest()
    actual_config = clustering_branch_config(plan, branch_id, config, source_sha256=source)
    labels = pd.Series([str(int(cell.split("_")[1]) % 3) for cell in spec["predicted_retained_ids"]],
                       index=spec["predicted_retained_ids"])
    digest = clustering_config_sha256(actual_config)
    return {"admitted": True, "run_status": "success", "source_sha256": source,
            "expected_source_sha256": source, "config_sha256": digest, "expected_config_sha256": digest,
            "scientific_config_sha256": clustering_scientific_config_sha256(actual_config),
            "raw_selected_ids": spec["selected_ids"], "selected_ids_sha256": spec["selected_ids_sha256"],
            "raw_counts_sha256": spec["raw_counts_sha256"], "feature_axis_sha256": spec["feature_axis_sha256"],
            "qc_parameters": spec["qc_parameters"], "clustering_parameters": spec["clustering_parameters"],
            "retained_ids": spec["predicted_retained_ids"], "clusters": labels}


def test_freezes_all_raw_source_masks_without_mutating_inputs_or_using_labels():
    payload, observations, config, plan = _fixture()
    before = deepcopy(config)
    validate_clustering_sensitivity_plan(plan)
    assert len(plan["branches"]) == 8
    assert [item["branch_id"] for item in plan["branches"][:4]] == [
        "resolution_0_5", "resolution_2_0", "clustering_seed_17", "uniform_sample_80pct"]
    sample = _spec(plan, "uniform_sample_80pct")
    expected = np.sort(np.random.default_rng(20261003).choice(80, size=64, replace=False))
    assert sample["selected_positions"] == expected.tolist()
    assert sample["selected_ids"] == payload.obs.index[expected].tolist()
    assert sample["raw_counts_sha256"] == raw_count_fingerprint(payload.layers["counts"][expected])
    assert sample["n_features"] == 30
    assert all(item["qc_parameters"] == plan["primary"]["qc_parameters"] for item in plan["branches"])
    assert len(sample["selected_ids"]) == 64  # Source sampling occurs before QC.
    assert plan["donor_perturbation"]["donor_independence_established"] is False
    assert config == before and observations == config["autonomous_analysis"]["observations"]


@pytest.mark.parametrize("donors, expected", [(0, "unavailable"), (1, "two to four"), (5, "two to four")])
def test_missing_or_out_of_bound_donors_are_not_invented_or_cherry_picked(donors, expected):
    _, _, _, plan = _fixture(donors=donors)
    assert len(plan["branches"]) == 4
    assert plan["donor_perturbation"]["status"] == "not_eligible"
    assert expected in plan["donor_perturbation"]["reason"]


def test_donor_leave_one_out_removes_complete_observed_identity_in_source_order():
    payload, _, _, plan = _fixture()
    for spec in plan["branches"][4:]:
        expected = payload.obs.index[payload.obs.donor_id != spec["omitted_donor"]].tolist()
        assert spec["selected_ids"] == expected
        assert spec["n_source_cells"] == 60
        assert spec["branch_id"] == "leave_out_" + spec["omitted_donor"]


def test_incomplete_donor_metadata_is_explicitly_ineligible():
    payload, _, config, _ = _fixture()
    payload.obs.loc["cell_7", "donor_id"] = None
    observation = observe_input(payload, question=config["autonomous_analysis"]["question"])
    config["autonomous_analysis"]["observations"] = observation
    config["capability_parameters"]["qc"] = {key: observation["qc"][key] for key in ("min_genes", "max_mito_pct")}
    plan = freeze_clustering_sensitivity_plan(payload=payload, observations=observation,
                                             source_sha256="a" * 64, baseline_config=config)
    assert plan["donor_perturbation"]["status"] == "not_eligible"
    assert "incomplete" in plan["donor_perturbation"]["reason"]
    assert len(plan["branches"]) == 4


def test_unsafe_donor_identifiers_receive_bounded_nonpath_slugs():
    payload, _, config, _ = _fixture(donors=2)
    payload.obs["donor_id"] = payload.obs["donor_id"].replace({"GSM2230757": "../escape/donor"})
    observation = observe_input(payload, question=config["autonomous_analysis"]["question"])
    config["autonomous_analysis"]["observations"] = observation
    plan = freeze_clustering_sensitivity_plan(payload=payload, observations=observation,
                                             source_sha256="a" * 64, baseline_config=config)
    spec = next(item for item in plan["branches"] if item["omitted_donor"] == "../escape/donor")
    assert spec["branch_id"].startswith("leave_out_donor_")
    assert "/" not in spec["branch_id"] and "\\" not in spec["branch_id"]


def test_case_distinct_donor_identities_do_not_collide_on_windows():
    payload, _, config, _ = _fixture(donors=2)
    payload.obs["donor_id"] = payload.obs["donor_id"].replace({"GSM2230757": "A", "GSM2230758": "a"})
    observation = observe_input(payload, question=config["autonomous_analysis"]["question"])
    config["autonomous_analysis"]["observations"] = observation
    plan = freeze_clustering_sensitivity_plan(payload=payload, observations=observation,
                                             source_sha256="a" * 64, baseline_config=config)
    ids = [item["branch_id"] for item in plan["branches"][4:]]
    assert len({name.casefold() for name in ids}) == 2
    assert all(name.startswith("leave_out_donor_") for name in ids)


def test_plan_rejects_nonexploratory_scope_and_unfrozen_primary_seed():
    payload, observations, config, _ = _fixture()
    for change in (lambda c: c["autonomous_analysis"].update(intent="condition_comparison"),
                   lambda c: c["capability_parameters"]["integration"].update(random_seed=17),
                   lambda c: c["capability_parameters"]["clustering"].update(resolution=.5),
                   lambda c: c["capability_parameters"]["clustering"].update(marker_reference={"ExpectedType": ["KnownGene"]}),
                   lambda c: c["analysis_extensions"].update(rna_velocity=True)):
        edited = deepcopy(config)
        change(edited)
        with pytest.raises(ValueError):
            freeze_clustering_sensitivity_plan(payload=payload, observations=observations,
                                                source_sha256="a" * 64, baseline_config=edited)


def test_not_estimable_source_subsample_is_recorded_without_rescue():
    _, _, config, plan = _fixture(donors=0, n_cells=3)
    spec = _spec(plan, "uniform_sample_80pct")
    assert spec["selected_positions"] == sorted(spec["selected_positions"])
    assert len(spec["selected_ids"]) == 2
    assert spec["planned_status"] == "not_estimable"
    assert spec["expected_observations"] is None
    with pytest.raises(ValueError, match="not-estimable"):
        clustering_branch_config(plan, "uniform_sample_80pct", config, source_sha256="b" * 64)


@pytest.mark.parametrize("source", ["X_integer_counts", "layers.counts"])
def test_subset_configuration_refreshes_canonical_source_facts_and_question_contract(source):
    payload, _, config, plan = _fixture(source=source)
    branch = clustering_branch_config(plan, "uniform_sample_80pct", config, source_sha256="b" * 64)
    settings = branch["autonomous_analysis"]
    assert settings["source_sha256"] == "b" * 64
    assert settings["observations"]["counts_source"] == "layers.counts"
    assert settings["observations"]["obs_column_map"]["donor"] == "donor"
    assert settings["requirement_contract"]["source_sha256"] == "b" * 64
    assert settings["observations"]["qc"] != config["autonomous_analysis"]["observations"]["qc"] or settings["observations"]["n_cells"] != payload.n_obs
    assert branch["capability_parameters"]["qc"] == config["capability_parameters"]["qc"]
    assert branch["capability_parameters"]["integration"] == config["capability_parameters"]["integration"]
    validate_clustering_branch_config(plan, "uniform_sample_80pct", branch)


@pytest.mark.parametrize("mutation", [
    lambda c: c["capability_parameters"]["qc"].update(min_genes=0),
    lambda c: c["capability_parameters"]["normalization"].update(n_top_genes=10),
    lambda c: c["capability_parameters"]["integration"].update(random_seed=17),
    lambda c: c["method_overrides"].update(integration="harmonypy_v1"),
    lambda c: c["analysis_extensions"].update(cluster_markers=False),
    lambda c: c["autonomous_analysis"]["observations"].update(n_cells=79),
    lambda c: c["autonomous_analysis"]["requirement_contract"].update(source_sha256="a" * 64),
])
def test_branch_validator_rejects_scientific_or_source_contract_drift(mutation):
    _, _, config, plan = _fixture()
    branch = clustering_branch_config(plan, "uniform_sample_80pct", config, source_sha256="b" * 64)
    mutation(branch)
    with pytest.raises(ValueError):
        validate_clustering_branch_config(plan, "uniform_sample_80pct", branch)


def test_no_self_hash_cycle_and_other_settings_remain_bound():
    _, _, config, plan = _fixture()
    digest = clustering_config_sha256(config)
    for name in ("clustering_sensitivity_plan", "clustering_sensitivity_branch", "qc_sensitivity_plan", "qc_sensitivity_branch"):
        config["autonomous_analysis"][name] = {"arbitrary": "self-reference"}
    assert clustering_config_sha256(config) == digest
    config["autonomous_analysis"]["other_setting"] = True
    assert clustering_config_sha256(config) != digest
    with pytest.raises(ValueError, match="Baseline configuration"):
        clustering_branch_config(plan, "resolution_0_5", config)


def test_plan_tamper_cannot_change_membership_or_parameter_delta():
    _, _, _, plan = _fixture()
    plan["branches"][0]["clustering_parameters"]["resolution"] = .8
    with pytest.raises(ValueError, match="plan identity"):
        validate_clustering_sensitivity_plan(plan)


def test_count_and_axis_hashes_detect_values_dtype_order_and_full_gene_universe():
    matrix = np.array([[2**24 + 1, 0, 2], [0, 4, 5]], dtype=np.int64)
    canonical = sparse.csr_matrix(matrix)
    assert raw_count_fingerprint(matrix) == raw_count_fingerprint(canonical)
    assert raw_count_fingerprint(matrix.astype(np.float32)) != raw_count_fingerprint(matrix)
    edited = matrix.copy()
    edited[0, 0] -= 1
    assert raw_count_fingerprint(edited) != raw_count_fingerprint(matrix)
    var = pd.DataFrame({"gene_name": ["A", "B", "C"]}, index=["id_a", "id_b", "id_c"])
    assert feature_axis_fingerprint(var, "var.gene_name") != feature_axis_fingerprint(var.iloc[::-1], "var.gene_name")
    assert feature_axis_fingerprint(var, "var.gene_name") != feature_axis_fingerprint(var.iloc[:2], "var.gene_name")


def test_label_permutations_and_series_reordering_do_not_change_paired_agreement():
    _, _, config, plan = _fixture()
    primary = _receipt(plan, config, "primary")
    branches = {item["branch_id"]: _receipt(plan, config, item["branch_id"]) for item in plan["branches"]}
    for receipt in branches.values():
        receipt["clusters"] = receipt["clusters"].map({"0": "C", "1": "A", "2": "B"}).iloc[::-1]
    compared = compare_clustering_sensitivity_runs(plan, baseline=primary, branches=branches)
    assert compared["status"] == "completed"
    assert all(item["shared_cell_ari"] == 1 and item["shared_cell_nmi"] == pytest.approx(1) for item in compared["branches"])
    assert compared["winner_selected"] is False and compared["primary_replaced"] is False
    assert compared["biological_stability_claim_supported"] is False


def test_ari_and_arithmetic_nmi_match_handworked_contingency_oracle():
    from eacbp.application.clustering_sensitivity import _paired_descriptions

    ids = list("abcdef")
    # Contingency [[2, 1], [0, 3]] gives pair sums 4, 6, 7,
    # chance expectation 42/15 and adjusted Rand index 12/37.
    primary = pd.Series(["A", "A", "A", "B", "B", "B"], index=ids)
    branch = pd.Series(["C", "C", "D", "D", "D", "D"], index=ids)
    result = _paired_descriptions(primary, branch, 20)
    expected_nmi = (3 * math.log(3) - 2 * math.log(2)) / (3 * math.log(3) + math.log(2))
    assert result["shared_cell_ari"] == pytest.approx(12 / 37)
    assert result["shared_cell_nmi"] == pytest.approx(expected_nmi)


def test_hungarian_matching_preserves_unmatched_clusters_and_small_support():
    _, _, config, plan = _fixture()
    primary = _receipt(plan, config, "primary")
    branch = _receipt(plan, config, "clustering_seed_17")
    primary["clusters"].iloc[:2] = "tiny"
    branch["clusters"][:] = "single"
    result = compare_clustering_sensitivity_runs(plan, baseline=primary, branches={"clustering_seed_17": branch})
    item = next(item for item in result["branches"] if item["branch_id"] == "clustering_seed_17")
    assert item["degenerate_partition"] is True
    assert len(item["cluster_matches"]) == 4
    assert sum(row["matching_status"] == "matched" for row in item["cluster_matches"]) == 1
    tiny = next(row for row in item["cluster_matches"] if row["primary_cluster"] == "tiny")
    assert tiny["small_support"] is True
    assert tiny["matching_status"] == "unmatched" and tiny["shared_jaccard"] == 0
    assert tiny["matched_branch_shared_cells"] == 0 and tiny["union_cells"] == 2
    assert item["primary_shared_fraction"] == 1 and item["branch_shared_fraction"] == 1


def test_two_single_cluster_partitions_have_explicit_degenerate_evidence_scope():
    _, _, config, plan = _fixture()
    primary = _receipt(plan, config, "primary")
    branch = _receipt(plan, config, "clustering_seed_17")
    primary["clusters"][:] = "primary_one"
    branch["clusters"][:] = "branch_one"
    result = compare_clustering_sensitivity_runs(plan, baseline=primary, branches={"clustering_seed_17": branch})
    item = next(item for item in result["branches"] if item["branch_id"] == "clustering_seed_17")
    assert item["shared_cell_ari"] == 1 and item["shared_cell_nmi"] == 1
    assert item["degenerate_partition"] is True
    assert "does not provide population-resolution evidence" in item["reason"]
    assert result["biological_stability_claim_supported"] is False


def test_unavoidable_zero_overlap_hungarian_pair_keeps_actual_partner_and_union():
    from eacbp.application.clustering_sensitivity import _paired_descriptions

    table = np.array([[3, 2, 0, 0], [2, 3, 0, 0], [2, 2, 0, 0], [0, 0, 3, 3]])
    main, branch = [], []
    for row in range(4):
        for column in range(4):
            main.extend([f"P{row}"] * int(table[row, column]))
            branch.extend([f"B{column}"] * int(table[row, column]))
    ids = [f"cell_{index}" for index in range(len(main))]
    result = _paired_descriptions(pd.Series(main, index=ids), pd.Series(branch, index=ids), 20)
    zeros = [row for row in result["cluster_matches"] if row["matching_status"] == "matched" and row["intersection_cells"] == 0]
    assert len(zeros) == 1
    row = zeros[0]
    assert row["matched_branch_cluster"] is not None and row["matched_branch_shared_cells"] > 0
    assert row["shared_jaccard"] == 0 and row["union_cells"] == row["shared_cells"] + row["matched_branch_shared_cells"]
    assert len(result["branch_cluster_matches"]) == 4
    assert all(item["matching_status"] == "matched" for item in result["branch_cluster_matches"])


def test_omitted_donor_only_cluster_remains_explicit_without_false_perfect_score():
    _, _, config, plan = _fixture()
    primary = _receipt(plan, config, "primary")
    branch_id = "leave_out_GSM2230757"
    branch = _receipt(plan, config, branch_id)
    omitted = set(plan["primary"]["predicted_retained_ids"]) - set(_spec(plan, branch_id)["predicted_retained_ids"])
    primary["clusters"].loc[list(omitted)] = "omitted_donor_cluster"
    result = compare_clustering_sensitivity_runs(plan, baseline=primary, branches={branch_id: branch})
    item = next(item for item in result["branches"] if item["branch_id"] == branch_id)
    row = next(row for row in item["cluster_matches"] if row["primary_cluster"] == "omitted_donor_cluster")
    assert row["shared_cells"] == 0 and row["shared_coverage"] == 0
    assert row["matching_status"] == "no_shared_cells" and row["shared_jaccard"] is None
    assert row["shared_recall"] is None and row["shared_precision"] is None
    assert row["matched_branch_shared_cells"] is None and row["union_cells"] is None


@pytest.mark.parametrize("mutation", [
    lambda r: r.update(raw_counts_sha256="f" * 64),
    lambda r: r.update(feature_axis_sha256="f" * 64),
    lambda r: r.update(selected_ids_sha256="f" * 64),
    lambda r: r.update(raw_selected_ids=r["raw_selected_ids"][::-1]),
    lambda r: r.update(retained_ids=r["retained_ids"][::-1]),
    lambda r: r.update(qc_parameters={"min_genes": 0, "max_mito_pct": 100.}),
    lambda r: r.update(clustering_parameters={"resolution": .9, "random_seed": 42, "n_neighbors": 15}),
    lambda r: r.update(scientific_config_sha256="f" * 64),
    lambda r: r.update(config_sha256="f" * 64),
    lambda r: r.update(source_sha256="f" * 64),
    lambda r: r.update(run_status="scientific_failure"),
    lambda r: r.update(clusters=pd.concat([r["clusters"], r["clusters"].iloc[:1]])),
    lambda r: r.update(clusters=r["clusters"].iloc[:-1]),
])
def test_admission_does_not_bless_changed_masks_counts_parameters_or_labels(mutation):
    _, _, config, plan = _fixture()
    primary = _receipt(plan, config, "primary")
    branch = _receipt(plan, config, "clustering_seed_17")
    mutation(branch)
    with pytest.raises(ValueError):
        compare_clustering_sensitivity_runs(plan, baseline=primary, branches={"clustering_seed_17": branch})


def test_unadmitted_primary_or_branch_does_not_supply_scientific_descriptions():
    _, _, config, plan = _fixture()
    primary = _receipt(plan, config, "primary")
    branch = _receipt(plan, config, "clustering_seed_17")
    branch["admitted"] = False
    result = compare_clustering_sensitivity_runs(plan, baseline=primary, branches={"clustering_seed_17": branch})
    item = next(item for item in result["branches"] if item["branch_id"] == "clustering_seed_17")
    assert item["status"] == "not_available" and "shared_cell_ari" not in item
    primary["admitted"] = False
    branch["admitted"] = True
    result = compare_clustering_sensitivity_runs(plan, baseline=primary, branches={"clustering_seed_17": branch})
    assert result["status"] == "partial"
    with pytest.raises(ValueError, match="Unexpected unfrozen"):
        compare_clustering_sensitivity_runs(plan, baseline=primary, branches={"optimised_winner": branch})
