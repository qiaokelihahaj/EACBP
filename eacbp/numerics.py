"""Small numerical helpers shared across EACBP analysis planes.

This module intentionally has no imports from capability or knowledge
packages.  Keeping the dependency direction neutral lets low level numerical
utilities be imported from any plane without triggering package initialisation
cycles.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np


def benjamini_hochberg(p_values: Sequence[float] | np.ndarray) -> np.ndarray:
    """Return BH adjusted p-values in the input order.

    Empty one-dimensional input returns an empty float array.  P-values are
    normalised explicitly before ranking for compatibility with historical
    callers: ``NaN`` and ``+inf`` are replaced by ``1``, ``-inf`` by ``0``,
    and finite values outside the mathematical ``[0, 1]`` interval are
    clipped to that interval.  These replacements are input-cleaning
    conventions and do not assign statistical meaning to invalid values.
    Ordinary finite p-values are unchanged, while adjusted values are always
    bounded by ``[0, 1]``.

    The helper intentionally accepts only a one-dimensional family.  This
    avoids silently flattening a matrix and accidentally changing the FDR
    family represented by a caller.
    """

    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1:
        raise ValueError("p_values must be a one-dimensional sequence")
    if values.size == 0:
        return np.empty(0, dtype=float)

    # Apply historical compatibility replacements before ranking.  The
    # replacements are numerical input cleaning, not interpretations of
    # invalid p-values as scientific evidence.
    cleaned = np.nan_to_num(values, nan=1.0, posinf=1.0, neginf=0.0)
    cleaned = np.clip(cleaned, 0.0, 1.0)
    order = np.argsort(cleaned, kind="stable")
    sorted_values = cleaned[order]
    adjusted_sorted = np.empty_like(sorted_values)
    running = 1.0
    count = len(sorted_values)
    for index in range(count - 1, -1, -1):
        running = min(running, sorted_values[index] * count / float(index + 1))
        adjusted_sorted[index] = min(1.0, max(0.0, running))

    result = np.empty_like(adjusted_sorted)
    result[order] = adjusted_sorted
    return result


__all__ = ["benjamini_hochberg"]
