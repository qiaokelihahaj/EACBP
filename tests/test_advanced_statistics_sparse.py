import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from unittest.mock import patch

from eacbp.auditor.advanced_statistics import _raw_counts, _reconstruct_input
from eacbp.capabilities.advanced_statistics import (
    AdvancedStatisticsInputError,
    _aggregate_pseudobulk,
    _prepare_pseudobulk,
    _validate_counts,
)
from eacbp.capabilities.sc_data import SCData


def _data(counts):
    metadata = []
    for donor in ("d1", "d2", "d3"):
        for condition in ("A", "B"):
            for cell in range(2):
                metadata.append({"donor_id": donor, "condition": condition, "cell": cell})
    return SCData(
        X=sparse.csr_matrix(np.zeros((len(metadata), 4), dtype=np.float32)),
        obs=pd.DataFrame(metadata),
        var=pd.DataFrame({"gene_name": ["g1", "g2", "g3", "g4"]}),
        layers={"counts": counts},
    )


def test_sparse_validation_and_audit_never_densify_cell_matrix():
    matrix = sparse.csr_matrix(np.arange(48, dtype=np.int64).reshape(12, 4))
    data = _data(matrix)
    with patch.object(sparse.csr_matrix, "toarray", side_effect=AssertionError("cell matrix was densified")):
        validated, _ = _validate_counts(data, {})
        audited = _raw_counts(data, {})
    assert sparse.issparse(validated)
    assert sparse.issparse(audited)
    assert validated.shape == audited.shape == (12, 4)


def test_full_independent_audit_reconstruction_never_densifies_cell_matrix():
    data = _data(sparse.csr_matrix(np.arange(48, dtype=np.int64).reshape(12, 4)))
    params = {"condition_a": "A", "condition_b": "B", "donor_col": "donor_id", "paired": True}
    with patch.object(sparse.csr_matrix, "toarray", side_effect=AssertionError("audit densified cell matrix")):
        reconstructed = _reconstruct_input(
            data, params, {"condition_a": "A", "condition_b": "B"}, "pydeseq2_pseudobulk_v1"
        )
    assert sparse.issparse(reconstructed.counts)
    assert reconstructed.design_rank == 4
    assert reconstructed.residual_df == 2


def test_sparse_pseudobulk_matches_dense_and_only_densifies_aggregated_rows():
    values = np.arange(1, 49, dtype=np.int64).reshape(12, 4)
    dense_data = _data(values)
    sparse_data = _data(sparse.csr_matrix(values))
    params = {"condition_a": "A", "condition_b": "B", "donor_col": "donor_id", "paired": True}
    dense = _prepare_pseudobulk(dense_data, params)

    called_shapes = []
    original_toarray = sparse.csr_matrix.toarray

    def tracked_toarray(matrix, *args, **kwargs):
        called_shapes.append(matrix.shape)
        assert matrix.shape != values.shape
        return original_toarray(matrix, *args, **kwargs)

    with patch.object(sparse.csr_matrix, "toarray", tracked_toarray):
        sparse_prepared = _prepare_pseudobulk(sparse_data, params)

    pd.testing.assert_frame_equal(sparse_prepared.counts, dense.counts)
    assert called_shapes == [(6, 4)]
    assert sparse_prepared.pseudobulk_dense_bytes == 6 * 4 * 8
    assert sparse_prepared.estimated_aggregation_working_bytes >= sparse_prepared.pseudobulk_dense_bytes


@pytest.mark.parametrize(
    "bad_values",
    [
        np.array([[1.5, 0.0]]),
        np.array([[-1, 0]], dtype=np.int64),
        np.array([[np.inf, 0.0]]),
        np.array([[float(2**63), 0.0]]),
        np.array([[2**63, 0]], dtype=np.uint64),
    ],
)
def test_sparse_count_validation_rejects_invalid_values(bad_values):
    data = SCData(
        X=np.zeros((1, 2), dtype=np.float32),
        obs=pd.DataFrame({"condition": ["A"], "donor": ["d1"]}),
        var=pd.DataFrame(index=["g1", "g2"]),
        layers={"counts": sparse.csr_matrix(bad_values)},
    )
    with pytest.raises(AdvancedStatisticsInputError):
        _validate_counts(data, {})


def test_sparse_duplicate_int8_counts_are_promoted_before_coalescing():
    duplicates = sparse.coo_matrix(
        (np.array([100, 100], dtype=np.int8), ([0, 0], [0, 0])), shape=(1, 1)
    )
    data = SCData(
        X=np.zeros((1, 1), dtype=np.float32),
        obs=pd.DataFrame({"condition": ["A"], "donor": ["d1"]}),
        var=pd.DataFrame(index=["g1"]),
        layers={"counts": duplicates},
    )
    validated, _ = _validate_counts(data, {})
    assert validated[0, 0] == 200


def test_pseudobulk_preflight_rejects_possible_int64_sum_overflow():
    obs = pd.DataFrame(
        {"donor": ["d1", "d1", "d2"], "condition": ["A", "A", "B"]}
    )
    counts = np.array([[np.iinfo(np.int64).max], [np.iinfo(np.int64).max], [0]], dtype=np.int64)
    with pytest.raises(AdvancedStatisticsInputError, match="may exceed int64"):
        _aggregate_pseudobulk(
            counts, obs, ["g1"], "condition", "donor", "A", "B", []
        )


def test_sparse_pseudobulk_resource_preflight_runs_before_densification():
    values = np.ones((12, 4), dtype=np.int64)
    data = _data(sparse.csr_matrix(values))
    params = {"condition_a": "A", "condition_b": "B", "donor_col": "donor_id", "paired": True}
    with patch("eacbp.capabilities.advanced_statistics._MAX_PSEUDOBULK_WORKING_BYTES", 1):
        with patch.object(sparse.csr_matrix, "toarray", side_effect=AssertionError("must preflight first")):
            with pytest.raises(AdvancedStatisticsInputError, match="resource limit"):
                _prepare_pseudobulk(data, params)
