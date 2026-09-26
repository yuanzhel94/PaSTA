"""On-demand Euclidean block-distance operations."""

from __future__ import annotations

import numpy as np

from ..types import BlockSetting


def blockij_euclidean_d2(coord_i, coord_j, coord_i_norm2=None, coord_j_norm2=None, is_diagonal_block: bool = False) -> np.ndarray:
    """Compute one squared Euclidean distance block without ``pdist``."""
    ci = np.asarray(coord_i, dtype=float)
    cj = np.asarray(coord_j, dtype=float)
    ni = np.sum(ci * ci, axis=1) if coord_i_norm2 is None else np.asarray(coord_i_norm2, dtype=float)
    nj = np.sum(cj * cj, axis=1) if coord_j_norm2 is None else np.asarray(coord_j_norm2, dtype=float)
    d2 = ni[:, None] + nj[None, :] - 2 * ci @ cj.T
    np.maximum(d2, 0, out=d2)
    if is_diagonal_block:
        np.fill_diagonal(d2, 0.0)
    return d2


def _block_distance(coord: np.ndarray, norm2: np.ndarray, blocks: BlockSetting, block_i: int, block_j: int) -> np.ndarray:
    si, ei = blocks.block_start[block_i], blocks.block_end[block_i]
    sj, ej = blocks.block_start[block_j], blocks.block_end[block_j]
    return np.sqrt(blockij_euclidean_d2(coord[si:ei], coord[sj:ej], norm2[si:ei], norm2[sj:ej], block_i == block_j))


def distance_limits_with_blocks_euclidean(coord: np.ndarray, coord_norm2: np.ndarray, setting) -> tuple[float, float]:
    """Compute global positive minimum and maximum Euclidean distance in blocks."""
    blocks = setting if isinstance(setting, BlockSetting) else __import__("pasta.config", fromlist=["set_blocks"]).set_blocks(setting.N, setting.block_size, setting.n_workers)
    dmin, dmax = np.inf, -np.inf
    for i, j in zip(blocks.block_i, blocks.block_j):
        d = _block_distance(coord, coord_norm2, blocks, int(i), int(j))
        values = d[np.triu_indices(d.shape[0], 1)] if i == j else d.ravel()
        values = values[values > 0]
        if values.size:
            dmin, dmax = min(dmin, float(values.min())), max(dmax, float(values.max()))
    if not np.isfinite(dmin):
        raise ValueError("coordinates contain no positive pairwise distances")
    return dmin, dmax


def within_range_count_euclidean(coord: np.ndarray, coord_norm2: np.ndarray, range_: float, setting) -> int:
    """Count ordered entries within an effective range blockwise."""
    blocks = setting if isinstance(setting, BlockSetting) else __import__("pasta.config", fromlist=["set_blocks"]).set_blocks(setting.N, setting.block_size, setting.n_workers)
    total = 0
    for i, j in zip(blocks.block_i, blocks.block_j):
        d = _block_distance(coord, coord_norm2, blocks, int(i), int(j))
        count = int(np.count_nonzero(d < range_))
        total += count if i == j else 2 * count
    return total
