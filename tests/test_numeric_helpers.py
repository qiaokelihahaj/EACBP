"""Focused regression tests for shared numerical helpers.

These tests intentionally exercise the small compatibility surface used by
the DEG, trajectory, spatial, and knowledge capabilities.  They do not test
the independent auditor implementations.
"""

import numpy as np
import subprocess
import sys
from pathlib import Path

from eacbp.capabilities.clustering import calculate_silhouette, simple_kmeans
from eacbp.capabilities.deg import benjamini_hochberg as deg_bh
from eacbp.capabilities.spatial import domain as spatial_domain
from eacbp.capabilities.spatial.autocorrelation import benjamini_hochberg as spatial_bh
from eacbp.capabilities.trajectory import _bh as trajectory_bh
from eacbp.knowledge.biological_db import BiologicalDBRetriever
from eacbp.numerics import benjamini_hochberg as shared_bh


def test_bh_known_values_are_adjusted_in_original_order():
    """The step-up correction is monotone after sorting and restores input order."""

    p_values = np.array([0.20, 0.01, 0.05], dtype=float)
    expected = np.array([0.20, 0.03, 0.075], dtype=float)

    assert deg_bh is shared_bh
    assert spatial_bh is shared_bh
    np.testing.assert_allclose(deg_bh(p_values), expected)
    np.testing.assert_allclose(spatial_bh(p_values), expected)
    np.testing.assert_allclose(trajectory_bh(p_values), expected)
    np.testing.assert_allclose(BiologicalDBRetriever.benjamini_hochberg(p_values.tolist()), expected)


def test_bh_empty_and_invalid_values_follow_compatibility_policy():
    """Empty input is an empty float array; invalid values are cleaned then clipped."""

    empty = deg_bh([])
    assert empty.shape == (0,)
    assert empty.dtype == np.dtype(float)

    invalid = np.array([np.nan, np.inf, -np.inf, -0.5, 1.5], dtype=float)
    # NaN/+inf -> 1 and -inf -> 0 retain their input positions.  Finite
    # out-of-range values are clipped before the BH ranking.
    expected = np.array([1.0, 1.0, 0.0, 0.0, 1.0], dtype=float)
    np.testing.assert_allclose(deg_bh(invalid), expected)
    assert np.all((deg_bh(invalid) >= 0.0) & (deg_bh(invalid) <= 1.0))


def test_bh_rejects_multidimensional_families_instead_of_flattening():
    """A matrix must not silently become one FDR family."""

    try:
        deg_bh(np.zeros((2, 2), dtype=float))
    except ValueError as exc:
        assert "one-dimensional" in str(exc)
    else:  # pragma: no cover - protects the policy if implementation changes
        raise AssertionError("multidimensional p-values were silently flattened")


def test_knowledge_imports_without_capabilities_package_cycle():
    """A fresh interpreter can import the knowledge database directly."""

    root = Path(__file__).resolve().parents[1]
    code = "from eacbp.knowledge.biological_db import BiologicalDBRetriever; print(BiologicalDBRetriever.__name__)"
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout.strip() == "BiologicalDBRetriever"


def test_kmeans_boundaries_are_deterministic_and_well_typed():
    """K-means clamps k to [1, n] and preserves deterministic integer labels."""

    X = np.array([[0.0], [1.0], [10.0]], dtype=np.float64)
    empty = simple_kmeans(np.empty((0, 1), dtype=np.float64), k=4)
    assert empty.shape == (0,)
    assert np.issubdtype(empty.dtype, np.integer)

    one = simple_kmeans(X, k=0, random_seed=7)
    assert one.dtype == np.dtype(int)
    assert np.array_equal(one, np.zeros(3, dtype=int))

    too_many = simple_kmeans(X, k=99, random_seed=7)
    assert too_many.shape == (3,)
    assert np.issubdtype(too_many.dtype, np.integer)
    assert np.array_equal(too_many, simple_kmeans(X, k=99, random_seed=7))
    assert set(too_many.tolist()).issubset({0, 1, 2})


def test_spatial_helpers_are_shared_with_float32_compatibility():
    """Spatial exports share K-means and retain silhouette's float32 contract."""

    assert spatial_domain.simple_kmeans is simple_kmeans

    X = np.array([[0.0, 0.0], [0.1, 0.0], [10.0, 10.0], [10.1, 10.0]], dtype=np.float64)
    labels = simple_kmeans(X, k=2, random_seed=42)
    score = spatial_domain.calculate_silhouette(X, labels, random_seed=42)
    expected = calculate_silhouette(X.astype(np.float32), labels, random_seed=42)
    assert isinstance(score, float)
    assert -1.0 <= score <= 1.0
    assert score == expected


def test_general_silhouette_keeps_float64_while_spatial_casts_float32():
    """The wrapper preserves the historical spatial precision boundary."""

    X = np.array([[1e8], [1e8 + 1], [1e8 + 128], [1e8 + 129]], dtype=np.float64)
    labels = np.array([0, 0, 1, 1], dtype=int)

    general = calculate_silhouette(X, labels)
    spatial = spatial_domain.calculate_silhouette(X, labels)

    assert spatial == calculate_silhouette(X.astype(np.float32), labels)
    assert spatial != general


def test_silhouette_single_cluster_and_empty_inputs_are_zero():
    """Degenerate label families have the historical neutral silhouette."""

    assert calculate_silhouette(np.empty((0, 2), dtype=np.float32), np.empty(0, dtype=int)) == 0.0
    assert calculate_silhouette(np.ones((3, 2), dtype=np.float32), np.zeros(3, dtype=int)) == 0.0
