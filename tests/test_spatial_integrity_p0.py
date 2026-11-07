"""Spatial conclusions require measurements, matched resources and section-local contacts."""

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.adapters.spacell import SpaCellAgentAdapter
from eacbp.capabilities.spatial import cci
from eacbp.capabilities.spatial.cci import CellCellCommunicationCapability, compute_spatial_cci
from eacbp.capabilities.spatial.autocorrelation import SpatialDEGCapability
from eacbp.capabilities.spatial.domain import (
    DENSE_GRAPH_LIMIT,
    SpatialDomainCapability,
    build_spatial_neighborhood_graph,
    resolve_spatial_sections,
    validate_spatial_graph,
)
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskStatus


def _data(*, coordinates=True, genes=None):
    coords = np.tile([[0., 0.], [1., 0.], [0., 1.]], (2, 1))
    return SCData(
        X=sparse.csr_matrix(np.array([[3, 1, 2], [2, 2, 1], [1, 4, 3]] * 2)),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B", "B", "B"],
                          "section_id": ["left"] * 3 + ["right"] * 3}),
        var=pd.DataFrame({"gene_name": genes or ["Apoe", "Trem2", "Gfap"]}),
        obsm={"spatial": coords} if coordinates else {},
    )


def _registered(tmp_path, data):
    registry = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    uri = "adata://spatial_integrity/input/v1"
    registry.register(uri_str=uri, payload=data.to_dict(), artifact_type=ArtifactType.SPATIAL_DATA,
                      study_id="spatial_integrity", created_by_task="input", operation="import")
    return registry, uri


def _contract(uri, capability, parameters=None):
    return TaskContract(task_id="spatial_check", capability=capability,
                        input_artifacts=[uri], parameters=parameters or {})


@pytest.mark.parametrize("capability", [SpatialDomainCapability, CellCellCommunicationCapability, SpatialDEGCapability])
def test_spatial_capabilities_reject_missing_measurements(tmp_path, capability):
    registry, uri = _registered(tmp_path, _data(coordinates=False))
    cap = capability()
    with pytest.raises(ValueError, match="Measured spatial coordinates not found"):
        cap.execute(_contract(uri, cap.capability_name), registry)
    assert not registry.exists("table://spatial_integrity/spatial_cci/v1")
    assert not registry.exists("adata://spatial_integrity/spatial_domains/v1")


def test_spacell_does_not_generate_grid_coordinates(tmp_path):
    registry, uri = _registered(tmp_path, _data(coordinates=False))
    adapter = SpaCellAgentAdapter()
    result = adapter.execute(_contract(uri, adapter.capability_name), registry)
    assert result.status == TaskStatus.EXECUTION_FAILURE
    assert "Measured spatial coordinates not found" in result.error_message
    assert result.output_artifacts == []


@pytest.mark.parametrize("pairs", [None, [], [("MissingLigand", "MissingReceptor", "resource")]])
def test_zero_lr_overlap_cannot_invent_gene_pairs(pairs):
    data = _data(genes=["gene_1", "gene_2", "gene_3"])
    sections, _ = resolve_spatial_sections(data)
    W, _, _ = build_spatial_neighborhood_graph(data.obsm["spatial"], section_labels=sections)
    with pytest.raises(ValueError, match="No ligand-receptor pairs.*not estimable"):
        compute_spatial_cci(data, W, lr_pairs=pairs, n_permutations=10)


def test_no_lr_match_registers_no_result(tmp_path):
    registry, uri = _registered(tmp_path, _data(genes=["gene_1", "gene_2", "gene_3"]))
    cap = CellCellCommunicationCapability()
    with pytest.raises(ValueError, match="No ligand-receptor pairs"):
        cap.execute(_contract(uri, cap.capability_name), registry)
    assert not registry.exists("table://spatial_integrity/spatial_cci/v1")


def test_identical_coordinates_in_two_sections_never_form_cross_edges():
    data = _data()
    sections, receipt = resolve_spatial_sections(data)
    graphs = build_spatial_neighborhood_graph(data.obsm["spatial"], k_neighbors=5, section_labels=sections)
    assert receipt["n_spatial_partitions"] == 2
    for graph in graphs:
        assert not np.any(graph[:3, 3:])
        assert not np.any(graph[3:, :3])
    assert np.allclose(graphs[0].sum(axis=1), 2)
    assert np.allclose(graphs[2].sum(axis=1), 1)


def test_large_partitioned_graph_remains_sparse():
    n = DENSE_GRAPH_LIMIT + 2
    half = n // 2
    coords = np.tile(np.column_stack([np.arange(half), np.zeros(half)]), (2, 1))
    sections = np.repeat(["section_1", "section_2"], half)
    W, D, W_norm = build_spatial_neighborhood_graph(coords, k_neighbors=3, section_labels=sections)
    for graph in (W, D, W_norm):
        assert sparse.isspmatrix_csr(graph)
        assert graph[:half, half:].nnz == 0
        assert graph[half:, :half].nnz == 0
        assert graph.nnz <= n * 6


@pytest.mark.parametrize("slot", ["obsp", "obsm", "uns"])
@pytest.mark.parametrize("capability", [CellCellCommunicationCapability, SpatialDEGCapability])
def test_supplied_cross_section_connectivity_is_rejected(tmp_path, slot, capability):
    data = _data()
    W = sparse.csr_matrix((np.ones(2), ([0, 3], [3, 0])), shape=(6, 6))
    getattr(data, slot)["spatial_connectivities"] = W
    registry, uri = _registered(tmp_path, data)
    cap = capability()
    with pytest.raises(ValueError, match="across observed sections"):
        cap.execute(_contract(uri, cap.capability_name), registry)


@pytest.mark.parametrize("sparse_graph", [False, True])
def test_graph_validation_rejects_bad_weights(sparse_graph):
    W = np.array([[0., np.nan], [1., 0.]])
    if sparse_graph:
        W = sparse.csr_matrix(W)
    with pytest.raises(ValueError, match="finite and nonnegative"):
        validate_spatial_graph(W, 2)


def test_permutations_preserve_each_section_label_distribution(monkeypatch):
    data = _data()
    sections, _ = resolve_spatial_sections(data)
    W, _, _ = build_spatial_neighborhood_graph(data.obsm["spatial"], section_labels=sections)
    original = cci.calculate_spatial_contact_density
    calls = []

    def capture(labels, unique_labels, graph):
        calls.append(np.asarray(labels).copy())
        return original(labels, unique_labels, graph)

    monkeypatch.setattr(cci, "calculate_spatial_contact_density", capture)
    result = compute_spatial_cci(data, W, n_permutations=20, random_seed=19)
    assert not result.empty
    assert len(calls) == 21
    for permuted in calls[1:]:
        for section in np.unique(sections):
            positions = np.flatnonzero(sections == section)
            assert sorted(permuted[positions]) == sorted(calls[0][positions])


def test_donor_local_section_ids_are_separate_coordinate_systems():
    data = _data()
    data.obs["section_id"] = "section_1"
    data.obs["donor_id"] = ["donor_a"] * 3 + ["donor_b"] * 3
    sections, receipt = resolve_spatial_sections(data)
    W, _, _ = build_spatial_neighborhood_graph(data.obsm["spatial"], section_labels=sections)
    assert receipt["n_spatial_partitions"] == 2
    assert not np.any(W[:3, 3:])


def test_multi_donor_input_cannot_assume_one_section_per_donor():
    data = _data()
    del data.obs["section_id"]
    data.obs["donor"] = ["a"] * 3 + ["b"] * 3
    with pytest.raises(ValueError, match="donor cannot be assumed to be one section"):
        resolve_spatial_sections(data)


def test_incomplete_section_metadata_is_not_pooled():
    data = _data()
    data.obs.loc[0, "section_id"] = None
    with pytest.raises(ValueError, match="complete and nonempty"):
        resolve_spatial_sections(data)


def test_observed_role_mapping_preserves_custom_section_identity():
    data = _data()
    data.obs.rename(columns={"section_id": "custom_slice"}, inplace=True)
    data.uns["autonomous_input_preparation"] = {"obs_column_map": {"section": "custom_slice"}}
    sections, receipt = resolve_spatial_sections(data)
    assert receipt["section_col"] == "custom_slice"
    assert len(np.unique(sections)) == 2


@pytest.mark.parametrize("capability", [SpatialDomainCapability, SpaCellAgentAdapter])
def test_domain_and_adapter_use_partitioned_measured_coordinates(tmp_path, capability):
    data = _data()
    data.obsm["spatial"] = np.column_stack([data.obsm["spatial"], np.arange(6)])
    registry, uri = _registered(tmp_path, data)
    cap = capability()
    result = cap.execute(_contract(uri, cap.capability_name, {"k_neighbors": 5, "n_domains": 2}), registry)
    assert result.status == TaskStatus.SUCCESS
    assert result.metrics["n_spatial_partitions"] == 2
    out_data = registry.load_scdata(result.output_artifacts[0])
    assert np.array_equal(out_data.obsm["spatial"], data.obsm["spatial"])
    W = out_data.obsp["spatial_connectivities"]
    assert np.count_nonzero(W[:3, 3:].toarray() if sparse.issparse(W) else W[:3, 3:]) == 0


def test_single_section_preserves_valid_cci_behavior(tmp_path):
    data = _data()
    del data.obs["section_id"]
    registry, uri = _registered(tmp_path, data)
    cap = CellCellCommunicationCapability()
    result = cap.execute(_contract(uri, cap.capability_name, {"n_permutations": 10}), registry)
    assert result.status == TaskStatus.SUCCESS
    assert result.metrics["n_spatial_partitions"] == 1
    assert result.metrics["permutation_scope"] == "within_observed_spatial_partition"
    _, table = registry.get(result.output_artifacts[0])
    assert set(table["ligand"]) == {"Apoe"}
    assert set(table["receptor"]) == {"Trem2"}


def test_multiple_sections_cannot_use_global_analytical_randomization(tmp_path):
    registry, uri = _registered(tmp_path, _data())
    cap = SpatialDEGCapability()
    with pytest.raises(ValueError, match="Multi-section.*section-aware inference"):
        cap.execute(_contract(uri, cap.capability_name), registry)
    assert not registry.exists("table://spatial_integrity/spatial_deg/v1")


def test_generic_replicate_metadata_does_not_become_biological_donor():
    data = _data()
    data.obs["replicate"] = ["library_a"] * 3 + ["library_b"] * 3
    _, receipt = resolve_spatial_sections(data)
    assert receipt["donor_col"] is None
    assert receipt["donor_role_source"] == "absent"


def test_mouse_alias_is_shared_with_biological_unit_resolution():
    data = _data()
    data.obs["section_id"] = "local_section"
    data.obs["mouse"] = ["mouse_a"] * 3 + ["mouse_b"] * 3
    sections, receipt = resolve_spatial_sections(data)
    assert receipt["donor_col"] == "mouse"
    assert receipt["n_spatial_partitions"] == 2
    assert len(np.unique(sections)) == 2


@pytest.mark.parametrize("mapping", ["explicit", "observed_role_map"])
def test_specimen_mapping_cannot_split_a_mouse_into_donors(mapping):
    data = _data()
    data.obs["mouse"] = "one_mouse"
    data.obs["sample"] = ["library_a"] * 3 + ["library_b"] * 3
    params = {"donor_col": "sample"} if mapping == "explicit" else {"obs_column_map": {"donor": "sample"}}
    with pytest.raises(ValueError, match="splits one biological individual"):
        resolve_spatial_sections(data, params)
