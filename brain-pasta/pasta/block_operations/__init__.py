"""Pair-block variogram and covariance operations."""

from __future__ import annotations

import numpy as np


def _accumulate_one(d, diff2, spec, is_diagonal_block):
    valid = (d > 0) & (d <= spec.variogram_dmax)
    if is_diagonal_block:
        valid &= np.triu(np.ones(d.shape, dtype=bool), 1)
    distance = d[valid]
    difference = diff2[valid]
    weight_sum = np.zeros(spec.M)
    weighted_sum = np.zeros(spec.M)
    if not distance.size:
        return weight_sum, weighted_sum
    nearest = np.rint((distance - spec.variogram_dmin) / spec.lag_sep).astype(int)
    for offset in spec.bin_idx_offset:
        bins = nearest + offset
        keep = (bins >= 0) & (bins < spec.M)
        if not np.any(keep):
            continue
        distances = distance[keep]
        bins = bins[keep]
        d2 = (distances - spec.h[bins]) ** 2
        inside = d2 <= spec.L2
        if not np.any(inside):
            continue
        bins, d2, values = bins[inside], d2[inside], difference[keep][inside]
        weights = np.exp(spec.weight_constant * d2)
        weight_sum += np.bincount(bins, weights=weights, minlength=spec.M)
        weighted_sum += np.bincount(bins, weights=weights * values, minlength=spec.M)
    return weight_sum, weighted_sum


def blockij_variogram_accumulation(d, xi, xj, yi, yj, spec, is_diagonal_block: bool):
    """Accumulate global variogram terms for both maps for one block pair."""
    dx = (np.asarray(xi)[:, None] - np.asarray(xj)[None, :]) ** 2
    dy = (np.asarray(yi)[:, None] - np.asarray(yj)[None, :]) ** 2
    weight, weighted_x = _accumulate_one(np.asarray(d), dx, spec, is_diagonal_block)
    _, weighted_y = _accumulate_one(np.asarray(d), dy, spec, is_diagonal_block)
    return weight, weighted_x, weighted_y


def blockij_variogram_accumulation_single_map(d, xi, xj, spec, is_diagonal_block: bool):
    """Accumulate one-map variogram terms for one block pair."""
    diff2 = (np.asarray(xi)[:, None] - np.asarray(xj)[None, :]) ** 2
    return _accumulate_one(np.asarray(d), diff2, spec, is_diagonal_block)


def blockij_stable_covariance(d, b, is_diagonal_block: bool):
    """Evaluate one stationary covariance block."""
    b = np.asarray(b, dtype=float)
    cov = b[0] * np.exp(-(np.asarray(d) / b[1]) ** b[2])
    if is_diagonal_block:
        np.fill_diagonal(cov, b[0] + b[3])
    return cov


def blockij_nonstationary_covariance(dij, bi, bj, dim: int, is_diagonal_block: bool):
    """Evaluate one PaSTA-NS process-convolution covariance block."""
    bi, bj = np.asarray(bi, dtype=float), np.asarray(bj, dtype=float)
    phii, phij = bi[:, 1], bj[:, 1]
    sig = (phii[:, None] ** 2 + phij[None, :] ** 2) / 2
    qij = np.asarray(dij) ** 2 / sig
    cov = ((phii[:, None] * phij[None, :]) / sig) ** (dim / 2) * np.sqrt(bi[:, 0, None] * bj[None, :, 0]) * np.exp(-(np.sqrt(qij) ** bi[0, 2]))
    if is_diagonal_block:
        np.fill_diagonal(cov, bi[:, 0] + bi[:, 3])
    return cov
