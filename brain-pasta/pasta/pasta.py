"""MATLAB-aligned, blockwise PaSTA entry point."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.stats import pearsonr

from .base import fit_stable_variogram, fit_stable_variogram_fixed_exponent, get_parc, nef2p, prepare_parcel_variogram_inputs, stats2nef, variogram_spec
from .block_operations import blockij_nonstationary_covariance, blockij_stable_covariance, blockij_variogram_accumulation, blockij_variogram_accumulation_single_map
from .config import pasta_setting, set_blocks
from .distance_euclidean import blockij_euclidean_d2
from .distance_variable import load_distance_blockij_in_memory, prepare_distance_in_memory
from .types import PaSTAFit, PaSTAModel


def _ordered_map(n_workers, fn, values):
    if n_workers > 1 and len(values) > 1:
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            return list(executor.map(fn, values))
    return [fn(value) for value in values]


def _global_variogram(x, y, get_distance, blocks, spec, n_workers):
    pairs = list(zip(blocks.block_i, blocks.block_j))
    def one(pair):
        i, j = map(int, pair)
        si, ei, sj, ej = blocks.block_start[i], blocks.block_end[i], blocks.block_start[j], blocks.block_end[j]
        return blockij_variogram_accumulation(get_distance(i, j), x[si:ei], x[sj:ej], y[si:ei], y[sj:ej], spec, i == j)
    weight, weighted_x, weighted_y = np.zeros(spec.M), np.zeros(spec.M), np.zeros(spec.M)
    for w, wx, wy in _ordered_map(n_workers, one, pairs):
        weight += w
        weighted_x += wx
        weighted_y += wy
    vx, vy = np.full(spec.M, np.nan), np.full(spec.M, np.nan)
    valid = weight > 0
    vx[valid], vy[valid] = 0.5 * weighted_x[valid] / weight[valid], 0.5 * weighted_y[valid] / weight[valid]
    return vx, vy


def _distance_limits(get_distance, blocks):
    dmin, dmax = np.inf, -np.inf
    for i, j in zip(blocks.block_i, blocks.block_j):
        d = get_distance(int(i), int(j))
        values = d[np.triu_indices(d.shape[0], 1)] if i == j else d.ravel()
        values = values[np.isfinite(values) & (values > 0)]
        if values.size:
            dmin, dmax = min(dmin, float(values.min())), max(dmax, float(values.max()))
    if not np.isfinite(dmin):
        raise ValueError("distances contain no finite positive pair")
    return dmin, dmax


def _within_range_count(range_, get_distance, blocks):
    count = 0
    for i, j in zip(blocks.block_i, blocks.block_j):
        n = int(np.count_nonzero(get_distance(int(i), int(j)) < range_))
        count += n if i == j else 2 * n
    return count


def _parcel_distance_limits(parc_idx, n_parc, get_distance, blocks):
    dmin, dmax = np.full(n_parc, np.inf), np.full(n_parc, -np.inf)
    for i, j in zip(blocks.block_i, blocks.block_j):
        i, j = int(i), int(j)
        si, ei, sj, ej = blocks.block_start[i], blocks.block_end[i], blocks.block_start[j], blocks.block_end[j]
        pi, pj, d = parc_idx[si:ei], parc_idx[sj:ej], get_distance(i, j)
        for parc in range(n_parc):
            values = d[np.ix_(pi == parc, pj == parc)].ravel()
            values = values[np.isfinite(values) & (values > 0)]
            if values.size:
                dmin[parc], dmax[parc] = min(dmin[parc], float(values.min())), max(dmax[parc], float(values.max()))
    if not np.all(np.isfinite(dmin)):
        raise ValueError("at least one parcel has no finite positive distances")
    return dmin, dmax


def _fit_parcel_models(x, parc_idx, n_parc, global_b, get_distance, blocks, setting):
    if n_parc == 1:
        return np.asarray(global_b, dtype=float).reshape(1, 4)
    dmin, dmax = _parcel_distance_limits(parc_idx, n_parc, get_distance, blocks)
    standardized, variances, specs = prepare_parcel_variogram_inputs(x, parc_idx, n_parc, dmin, dmax, setting)
    weights, weighted = [np.zeros(s.M) for s in specs], [np.zeros(s.M) for s in specs]
    for i, j in zip(blocks.block_i, blocks.block_j):
        i, j = int(i), int(j)
        si, ei, sj, ej = blocks.block_start[i], blocks.block_end[i], blocks.block_start[j], blocks.block_end[j]
        pi, pj, d = parc_idx[si:ei], parc_idx[sj:ej], get_distance(i, j)
        for parc, spec in enumerate(specs):
            mi, mj = pi == parc, pj == parc
            if mi.any() and mj.any():
                w, wd = blockij_variogram_accumulation_single_map(d[np.ix_(mi, mj)], standardized[si:ei][mi], standardized[sj:ej][mj], spec, i == j)
                weights[parc] += w
                weighted[parc] += wd
    bs = np.empty((n_parc, 4))
    for parc, spec in enumerate(specs):
        v = np.full(spec.M, np.nan)
        valid = weights[parc] > 0
        v[valid] = 0.5 * weighted[parc][valid] / weights[parc][valid]
        b, _ = fit_stable_variogram_fixed_exponent(v, spec.h, global_b[2], setting.nugget)
        b[0], b[3] = b[0] * variances[parc], b[3] * variances[parc]
        bs[parc] = b
    return bs


def _prepare_stats(get_distance, blocks, bx, by, parc_x, parc_y, dim, n_workers, nonstationary):
    trace_x = float(np.sum(bx[parc_x, 0] + bx[parc_x, 3]))
    trace_y = float(np.sum(by[parc_y, 0] + by[parc_y, 3]))
    pairs = list(zip(blocks.block_i, blocks.block_j))
    def one(pair):
        i, j = map(int, pair)
        si, ei, sj, ej = blocks.block_start[i], blocks.block_end[i], blocks.block_start[j], blocks.block_end[j]
        d = get_distance(i, j)
        if nonstationary:
            cx = blockij_nonstationary_covariance(d, bx[parc_x[si:ei]], bx[parc_x[sj:ej]], dim, i == j)
            cy = blockij_nonstationary_covariance(d, by[parc_y[si:ei]], by[parc_y[sj:ej]], dim, i == j)
        else:
            cx, cy = blockij_stable_covariance(d, bx[0], i == j), blockij_stable_covariance(d, by[0], i == j)
        return i, j, np.sum(cx, 1), np.sum(cy, 1), None if i == j else np.sum(cx, 0), None if i == j else np.sum(cy, 0), float(np.dot(cx.ravel(), cy.ravel())) * (1 if i == j else 2)
    rowsum_x, rowsum_y, product = np.zeros(blocks.N), np.zeros(blocks.N), 0.0
    for i, j, rx, ry, cx, cy, value in _ordered_map(n_workers, one, pairs):
        si, ei = blocks.block_start[i], blocks.block_end[i]
        rowsum_x[si:ei] += rx
        rowsum_y[si:ei] += ry
        if i != j:
            sj, ej = blocks.block_start[j], blocks.block_end[j]
            rowsum_x[sj:ej] += cx
            rowsum_y[sj:ej] += cy
        product += value
    return rowsum_x, rowsum_y, trace_x, trace_y, product


def _filtered_parc(parc, valid, N, name):
    if parc is None or isinstance(parc, str):
        return parc
    value = np.asarray(parc)
    if value.ndim != 1 or value.size != N:
        raise ValueError(f"{name} must have one entry per original observation")
    return value[valid]


def pasta(x, y, coord, *, M=None, M_scale=3, qd=0.7, nugget=True, kernel_scale=4.0, block_size=2000, n_workers=1, xparc=None, yparc=None, max_clusters=10, min_clusters=1, min_cluster_size=500, D=None, dim=None, random_state=None) -> PaSTAFit:
    """Run memory-efficient PaSTA or PaSTA-NS.

    The keyword interface and returned fields match ``matlab_PaSTA/pasta.m``.
    ``D`` accepts a full matrix or a NumPy-order packed strict upper triangle.
    """
    x, y, coordinate = np.asarray(x, dtype=float).reshape(-1), np.asarray(y, dtype=float).reshape(-1), np.asarray(coord, dtype=float)
    if x.size != y.size:
        raise ValueError("x and y must have the same length")
    N_in = x.size
    if coordinate.ndim != 2 or coordinate.shape[0] != N_in:
        raise ValueError("coord must be a real numeric matrix with one row per observation")
    if random_state is not None and (not isinstance(random_state, (int, np.integer)) or random_state < 0):
        raise ValueError("random_state must be a non-negative integer or None")
    valid = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(coordinate), axis=1)
    xparc, yparc = _filtered_parc(xparc, valid, N_in, "xparc"), _filtered_parc(yparc, valid, N_in, "yparc")
    x, y, coordinate = x[valid], y[valid], coordinate[valid]
    if x.size < 2:
        raise ValueError("fewer than two finite observations remain")
    std_x, std_y = float(np.std(x, ddof=1)), float(np.std(y, ddof=1))
    if not np.isfinite(std_x) or std_x <= 0:
        raise ValueError("x has zero or invalid variance")
    if not np.isfinite(std_y) or std_y <= 0:
        raise ValueError("y has zero or invalid variance")
    rX, p_naive = pearsonr(x, y)
    x, y = (x - x.mean()) / std_x, (y - y.mean()) / std_y
    dim = coordinate.shape[1] if dim is None else dim
    if not isinstance(dim, (int, np.integer)) or dim < 1:
        raise ValueError("dim must be a positive integer")
    setting = pasta_setting(x.size, M=M, M_scale=M_scale, qd=qd, nugget=nugget, kernel_scale=kernel_scale, block_size=block_size, n_workers=n_workers, xparc=xparc, yparc=yparc, max_clusters=max_clusters, min_clusters=min_clusters, min_cluster_size=min_cluster_size)
    blocks = set_blocks(setting.N, setting.block_size, setting.n_workers)
    if D is None:
        norm2 = np.sum(coordinate * coordinate, axis=1)
        def get_distance(i, j):
            si, ei, sj, ej = blocks.block_start[i], blocks.block_end[i], blocks.block_start[j], blocks.block_end[j]
            return np.sqrt(blockij_euclidean_d2(coordinate[si:ei], coordinate[sj:ej], norm2[si:ei], norm2[sj:ej], i == j))
    else:
        memory = prepare_distance_in_memory(D, N_in, np.flatnonzero(valid))
        def get_distance(i, j):
            return load_distance_blockij_in_memory(memory, i, j, blocks)
    dmin, dmax = _distance_limits(get_distance, blocks)
    spec = variogram_spec(dmin, dmax, setting.M, setting.qd, setting.kernel_scale, setting.nugget)
    vx, vy = _global_variogram(x, y, get_distance, blocks, spec, setting.n_workers)
    global_bx, fx = fit_stable_variogram(vx, spec.h, setting.nugget)
    global_by, fy = fit_stable_variogram(vy, spec.h, setting.nugget)
    range_x = global_bx[1] * 2.996 ** (1 / global_bx[2]) if xparc is not None else np.nan
    range_y = global_by[1] * 2.996 ** (1 / global_by[2]) if yparc is not None else np.nan
    rng = np.random.RandomState(random_state) if random_state is not None else None
    x_seed = None if rng is None else int(rng.randint(0, np.iinfo(np.int32).max))
    y_seed = None if rng is None else int(rng.randint(0, np.iinfo(np.int32).max))
    parc_x, _, n_parc_x = get_parc(xparc, coordinate, _within_range_count(range_x, get_distance, blocks) / setting.N - 1 if xparc is not None else np.nan, setting, 1, x_seed)
    parc_y, _, n_parc_y = get_parc(yparc, coordinate, _within_range_count(range_y, get_distance, blocks) / setting.N - 1 if yparc is not None else np.nan, setting, 2, y_seed)
    bs_x = _fit_parcel_models(x, parc_x, n_parc_x, global_bx, get_distance, blocks, setting)
    bs_y = _fit_parcel_models(y, parc_y, n_parc_y, global_by, get_distance, blocks, setting)
    rowsum_x, rowsum_y, trace_x, trace_y, cov_product = _prepare_stats(get_distance, blocks, bs_x, bs_y, parc_x, parc_y, int(dim), setting.n_workers, n_parc_x > 1 or n_parc_y > 1)
    nef, run_status = stats2nef(rowsum_x, rowsum_y, trace_x, trace_y, cov_product)
    return PaSTAFit(nef2p(rX, nef) if run_status else np.nan, float(rX), float(nef), run_status, np.array([n_parc_x, n_parc_y]), float(p_naive), PaSTAModel(global_bx, bs_x, parc_x, n_parc_x, vx, spec.h, fx), PaSTAModel(global_by, bs_y, parc_y, n_parc_y, vy, spec.h, fy))
