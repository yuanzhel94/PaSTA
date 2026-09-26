"""Configuration helpers corresponding to ``matlab_PaSTA/config``."""

from __future__ import annotations

import numpy as np

from ..types import BlockSetting, PaSTASetting


def set_blocks(N: int, block_size: int | None = None, n_workers: int = 1) -> BlockSetting:
    """Create zero-based upper-triangular block-pair indices."""
    if not isinstance(N, (int, np.integer)) or N < 1:
        raise ValueError("N must be a positive integer")
    if block_size is None:
        block_size = N
    if not isinstance(block_size, (int, np.integer)) or block_size < 1:
        raise ValueError("block_size must be a positive integer or None")
    if not isinstance(n_workers, (int, np.integer)) or n_workers < 1:
        raise ValueError("n_workers must be a positive integer")
    N = int(N)
    block_size = min(int(block_size), N)
    starts = np.arange(0, N, block_size, dtype=int)
    ends = np.minimum(starts + block_size, N)
    pairs = [(i, j) for i in range(starts.size) for j in range(i, starts.size)]
    block_i = np.fromiter((p[0] for p in pairs), dtype=int)
    block_j = np.fromiter((p[1] for p in pairs), dtype=int)
    return BlockSetting(N, block_size, int(n_workers), starts, ends, block_i, block_j, block_i == block_j)


def pasta_setting(
    N: int,
    *,
    M: int | None = None,
    M_scale: int = 3,
    qd: float = 0.7,
    nugget: bool = True,
    kernel_scale: float = 4.0,
    block_size: int | None = 2000,
    n_workers: int = 1,
    xparc=None,
    yparc=None,
    max_clusters: int = 10,
    min_clusters: int = 1,
    min_cluster_size: int = 500,
) -> PaSTASetting:
    """Validate PaSTA numerical settings and set MATLAB-equivalent defaults."""
    if not isinstance(N, (int, np.integer)) or N < 2:
        raise ValueError("at least two valid observations are required")
    for name, value in (("M_scale", M_scale), ("max_clusters", max_clusters), ("min_clusters", min_clusters), ("min_cluster_size", min_cluster_size), ("n_workers", n_workers)):
        if not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if M is not None and (not isinstance(M, (int, np.integer)) or M < 2):
        raise ValueError("M must be an integer of at least 2 or None")
    if not np.isfinite(qd) or not 0 < qd <= 1:
        raise ValueError("qd must be finite and in (0, 1]")
    if not np.isfinite(kernel_scale) or kernel_scale <= 0:
        raise ValueError("kernel_scale must be a positive finite number")
    if min_clusters > max_clusters:
        raise ValueError("min_clusters cannot exceed max_clusters")
    N = int(N)
    M_scale = int(M_scale)
    return PaSTASetting(
        N=N,
        M=int(M) if M is not None else M_scale * int(np.ceil(np.sqrt(N))),
        M_scale=M_scale,
        qd=float(qd),
        nugget=bool(nugget),
        kernel_scale=float(kernel_scale),
        block_size=N if block_size is None else min(int(block_size), N),
        n_workers=int(n_workers),
        xparc=xparc,
        yparc=yparc,
        max_clusters=int(max_clusters),
        min_clusters=int(min_clusters),
        min_cluster_size=int(min_cluster_size),
    )
