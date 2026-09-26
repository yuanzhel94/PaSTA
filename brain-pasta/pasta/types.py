"""Public data structures used by the MATLAB-aligned PaSTA API."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal

import numpy as np


@dataclass(frozen=True)
class BlockSetting:
    """Block partition used to avoid dense pairwise allocations."""

    N: int
    block_size: int
    n_workers: int
    block_start: np.ndarray
    block_end: np.ndarray
    block_i: np.ndarray
    block_j: np.ndarray
    is_diagonal_block: np.ndarray

    @property
    def n_blocks(self) -> int:
        return int(self.block_start.size)

    @property
    def n_unique_block_pairs(self) -> int:
        return int(self.block_i.size)

    @property
    def parallel(self) -> bool:
        return self.n_workers > 1 and self.n_unique_block_pairs > 1


@dataclass(frozen=True)
class PaSTASetting:
    """Validated numerical settings corresponding to MATLAB ``setting``."""

    N: int
    M: int
    M_scale: int
    qd: float
    nugget: bool
    kernel_scale: float
    block_size: int
    n_workers: int
    xparc: object
    yparc: object
    max_clusters: int
    min_clusters: int
    min_cluster_size: int


@dataclass(frozen=True)
class VariogramSpec:
    """Lag-grid and truncated Gaussian-kernel specification."""

    M: int
    qd: float
    kernel_scale: float
    nugget: bool
    dmin: float
    dmax: float
    variogram_dmin: float
    variogram_dmax: float
    h: np.ndarray
    lag_sep: float
    delta: float
    sigma: float
    L: float
    L2: float
    weight_constant: float
    kernel_nh: int
    bin_idx_offset: np.ndarray


@dataclass(frozen=True)
class DistanceMemory:
    """An in-memory full or NumPy-packed strict-upper-triangle distance input."""

    data: np.ndarray
    storage: Literal["matrix", "upper"]
    N: int
    indices: np.ndarray | None = None


@dataclass(frozen=True)
class StableVariogram:
    """Callable stable semivariogram, replacing MATLAB's function handle."""

    b: np.ndarray

    def __call__(self, distance: np.ndarray | float) -> np.ndarray | float:
        sill, scale, exponent, nugget = self.b
        return sill * (1.0 - np.exp(-(np.asarray(distance) / scale) ** exponent)) + nugget


@dataclass(frozen=True)
class PaSTAModel:
    """Per-map model fields corresponding to MATLAB ``model_x`` / ``model_y``."""

    global_b: np.ndarray
    parcel_b: np.ndarray
    parc_idx: np.ndarray
    n_parc: int
    v: np.ndarray
    h: np.ndarray
    f: Callable[[np.ndarray | float], np.ndarray | float]


@dataclass(frozen=True)
class PaSTAFit:
    """Result of :func:`pasta.pasta`, equivalent to MATLAB ``pasta_fit``."""

    pef: float
    rX: float
    nef: float
    run_status: bool
    n_parc: np.ndarray
    p_naive: float
    model_x: PaSTAModel
    model_y: PaSTAModel
