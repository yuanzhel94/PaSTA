"""Blockwise access to in-memory non-Euclidean distance inputs."""

from __future__ import annotations

import numpy as np

from ..types import BlockSetting, DistanceMemory


def prepare_distance_in_memory(distance_data, N: int, indices: np.ndarray | None = None) -> DistanceMemory:
    """Validate a full matrix or NumPy-packed strict-upper-triangle vector."""
    data = np.asarray(distance_data)
    if not np.issubdtype(data.dtype, np.number) or np.iscomplexobj(data):
        raise ValueError("D must be a real numeric array")
    if data.ndim == 2 and data.shape == (N, N):
        return DistanceMemory(data, "matrix", int(N), indices)
    if data.ndim == 1 and data.size == N * (N - 1) // 2:
        return DistanceMemory(data, "upper", int(N), indices)
    raise ValueError("D must have shape (N, N) or be a NumPy-order strict upper-triangle vector")


def _upper_index(N: int, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Return NumPy ``triu_indices`` packed locations for ``i < j``."""
    return i * (2 * N - i - 1) // 2 + (j - i - 1)


def load_distance_blockij_in_memory(distance_memory: DistanceMemory, block_i: int, block_j: int, distance_blocks: BlockSetting) -> np.ndarray:
    """Extract one complete symmetric diagonal or upper off-diagonal distance block."""
    if block_j < block_i:
        raise ValueError("block_j must be greater than or equal to block_i")
    rs, re = distance_blocks.block_start[block_i], distance_blocks.block_end[block_i]
    cs, ce = distance_blocks.block_start[block_j], distance_blocks.block_end[block_j]
    rows = np.arange(rs, re, dtype=np.int64)
    cols = np.arange(cs, ce, dtype=np.int64)
    if distance_memory.indices is not None:
        rows, cols = distance_memory.indices[rows], distance_memory.indices[cols]
    if distance_memory.storage == "matrix":
        return np.asarray(distance_memory.data[np.ix_(rows, cols)], dtype=float)
    output = np.zeros((rows.size, cols.size), dtype=float)
    if block_i == block_j:
        ii, jj = np.triu_indices(rows.size, k=1)
        gi, gj = rows[ii], rows[jj]
        values = np.asarray(distance_memory.data[_upper_index(distance_memory.N, gi, gj)], dtype=float)
        output[ii, jj] = values
        output[jj, ii] = values
    else:
        gi = rows[:, None]
        gj = cols[None, :]
        output[:] = np.asarray(distance_memory.data[_upper_index(distance_memory.N, gi, gj)], dtype=float)
    return output


def distance_limits_with_blocks_in_memory(distance_memory: DistanceMemory, blocks: BlockSetting):
    """Compute positive minimum and maximum distance without dense expansion."""
    dmin, dmax = np.inf, -np.inf
    for i, j in zip(blocks.block_i, blocks.block_j):
        d = load_distance_blockij_in_memory(distance_memory, int(i), int(j), blocks)
        if i == j:
            values = d[np.triu_indices(d.shape[0], k=1)]
        else:
            values = d.ravel()
        finite = values[np.isfinite(values) & (values > 0)]
        if finite.size:
            dmin = min(dmin, float(finite.min()))
            dmax = max(dmax, float(finite.max()))
    if not np.isfinite(dmin) or not np.isfinite(dmax):
        raise ValueError("D contains no finite positive distances")
    return dmin, dmax


def within_range_count_in_memory(ranges: np.ndarray, distance_memory: DistanceMemory, blocks: BlockSetting) -> np.ndarray:
    """Count ordered entries with distance below each supplied effective range."""
    ranges = np.asarray(ranges, dtype=float).reshape(-1)
    counts = np.zeros(ranges.size, dtype=np.int64)
    for i, j in zip(blocks.block_i, blocks.block_j):
        d = load_distance_blockij_in_memory(distance_memory, int(i), int(j), blocks)
        for k, value in enumerate(ranges):
            count = int(np.count_nonzero(d < value))
            if i != j:
                count *= 2
            counts[k] += count
    return counts
