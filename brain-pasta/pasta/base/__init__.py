"""Core model helpers corresponding to ``matlab_PaSTA/base``."""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import t as t_dist
from sklearn.cluster import KMeans

from ..types import PaSTASetting, StableVariogram, VariogramSpec


def variogram_spec(dmin: float, dmax: float, M: int, qd: float, kernel_scale: float, nugget: bool) -> VariogramSpec:
    """Create the lag grid and truncated Gaussian smoothing kernel."""
    dmin, dmax = float(dmin), float(dmax)
    vmax = qd * dmax
    if not np.isfinite(dmin) or not np.isfinite(dmax) or dmin <= 0 or vmax <= dmin:
        raise ValueError("qd*dmax must be greater than a positive finite dmin")
    h = np.linspace(dmin, vmax, int(M))
    lag_sep = (vmax - dmin) / (int(M) - 1)
    delta = lag_sep / 2
    sigma = 6 * delta
    L = float(kernel_scale) * sigma / 2.68
    kernel_nh = int(np.ceil(L / lag_sep))
    return VariogramSpec(int(M), float(qd), float(kernel_scale), bool(nugget), dmin, dmax, dmin, vmax, h, lag_sep, delta, sigma, L, L * L, -2.68**2 / (2 * sigma**2), kernel_nh, np.arange(-kernel_nh, kernel_nh + 1, dtype=int))


def _stable(distance, sill, scale, exponent, nugget):
    return sill * (1.0 - np.exp(-(distance / scale) ** exponent)) + nugget


def fit_stable_variogram(v: np.ndarray, h: np.ndarray, nugget: bool = True):
    """Fit MATLAB's stable variogram and return ``(b, f)``."""
    v, h = np.asarray(v, dtype=float), np.asarray(h, dtype=float)
    valid = np.isfinite(h) & np.isfinite(v) & (h > 0) & (v >= 0)
    if valid.sum() < 2:
        raise ValueError("variogram has fewer than two valid lag bins")
    h_fit, v_fit = h[valid], v[valid]
    v_max, h_min = float(np.max(v_fit)), float(np.min(h_fit))
    if v_max <= 0:
        raise ValueError("variogram has zero variance")
    if nugget:
        b, _ = curve_fit(_stable, h_fit, v_fit, p0=(v_max, h_min, 1.0, 0.0), bounds=((0, 0, 0, 0), (2 * v_max, np.inf, 2, 0.5 * v_max)))
    else:
        def f_no_nugget(d, sill, scale, exponent):
            return _stable(d, sill, scale, exponent, 0.0)
        p, _ = curve_fit(f_no_nugget, h_fit, v_fit, p0=(v_max, h_min, 1.0), bounds=((0, 0, 0), (2 * v_max, np.inf, 2)))
        b = np.array([*p, 0.0])
    b = np.asarray(b, dtype=float)
    return b, StableVariogram(b)


def fit_stable_variogram_fixed_exponent(v: np.ndarray, h: np.ndarray, exponent: float, nugget: bool = True):
    """Fit the stable model while holding its exponent fixed."""
    v, h = np.asarray(v, dtype=float), np.asarray(h, dtype=float)
    valid = np.isfinite(h) & np.isfinite(v) & (h > 0) & (v >= 0)
    if valid.sum() < 2:
        raise ValueError("variogram has fewer than two valid lag bins")
    h_fit, v_fit = h[valid], v[valid]
    v_max, h_min = float(np.max(v_fit)), float(np.min(h_fit))
    if v_max <= 0:
        raise ValueError("variogram has zero variance")
    if nugget:
        def fixed(d, sill, scale, nug):
            return _stable(d, sill, scale, exponent, nug)
        p, _ = curve_fit(fixed, h_fit, v_fit, p0=(v_max, h_min, 0.0), bounds=((0, 0, 0), (2 * v_max, np.inf, 0.5 * v_max)))
        b = np.array([p[0], p[1], exponent, p[2]])
    else:
        def fixed(d, sill, scale):
            return _stable(d, sill, scale, exponent, 0.0)
        p, _ = curve_fit(fixed, h_fit, v_fit, p0=(v_max, h_min), bounds=((0, 0), (2 * v_max, np.inf)))
        b = np.array([p[0], p[1], exponent, 0.0])
    return b, StableVariogram(b)


def get_parc(parc_arg, coord: np.ndarray, n_points_from_range: float, setting: PaSTASetting, map_idx: int, random_state: int | None = None):
    """Return zero-based parcel indices, labels, and parcel count."""
    if parc_arg is None:
        return np.zeros(setting.N, dtype=int), np.empty(0, dtype=int), 1
    n_points = max(float(n_points_from_range), setting.min_cluster_size)
    n_parc_auto = max(min(int(np.floor(setting.N / n_points)), setting.max_clusters), setting.min_clusters)
    if isinstance(parc_arg, str):
        if parc_arg != "auto":
            raise ValueError(f"invalid parcellation input for map {map_idx}")
        if n_parc_auto == 1:
            return np.zeros(setting.N, dtype=int), np.array([0]), 1
        labels = KMeans(n_clusters=n_parc_auto, init="k-means++", n_init=1, random_state=random_state).fit_predict(coord)
        unique, inverse = np.unique(labels, return_inverse=True)
        return inverse.astype(int), unique, int(unique.size)
    labels = np.asarray(parc_arg)
    if labels.ndim != 1 or labels.size != setting.N:
        raise ValueError(f"user specified parcellation for map {map_idx} must have one entry per valid observation")
    unique, inverse = np.unique(labels, return_inverse=True)
    if unique.size > n_parc_auto:
        warnings.warn(f"data No.{map_idx}: specified number of parcs {unique.size} exceeds data-derived maximum {n_parc_auto}", UserWarning, stacklevel=2)
    return inverse.astype(int), unique, int(unique.size)


def prepare_parcel_variogram_inputs(x: np.ndarray, parc_idx: np.ndarray, n_parc: int, dmin: np.ndarray, dmax: np.ndarray, setting: PaSTASetting):
    """Standardize each parcel and prepare its local variogram specification."""
    standardized = np.empty_like(x, dtype=float)
    variances = np.empty(n_parc, dtype=float)
    specs = []
    for parc in range(n_parc):
        selected = parc_idx == parc
        values = x[selected]
        variance = float(np.var(values, ddof=1))
        if not np.isfinite(variance) or variance <= 0:
            raise ValueError(f"Parcel {parc} has zero or invalid variance")
        standardized[selected] = (values - values.mean()) / np.sqrt(variance)
        variances[parc] = variance
        specs.append(variogram_spec(dmin[parc], dmax[parc], setting.M_scale * int(np.ceil(np.sqrt(values.size))), setting.qd, setting.kernel_scale, setting.nugget))
    return standardized, variances, specs


def stats2nef(rowsum_x: np.ndarray, rowsum_y: np.ndarray, trace_x: float, trace_y: float, cov_prodsum: float):
    """Compute effective sample size from blockwise covariance statistics."""
    N = rowsum_x.size
    sum_x, sum_y = float(np.sum(rowsum_x)), float(np.sum(rowsum_y))
    centered_trace_x = trace_x - sum_x / N
    centered_trace_y = trace_y - sum_y / N
    cov_prod = cov_prodsum - 2 * float(np.dot(rowsum_x, rowsum_y)) / N + (sum_x * sum_y) / N**2
    nef = float(np.real(centered_trace_x * centered_trace_y / cov_prod + 1)) if cov_prod != 0 else np.nan
    return nef, bool(np.isfinite(nef) and np.isfinite(cov_prod) and cov_prod > 0 and nef > 2)


def nef2p(rX: float, nef: float) -> float:
    """Two-sided correlation p-value using PaSTA effective sample size."""
    df = max(0.0, float(nef) - 2.0)
    if df == 0 or not np.isfinite(df):
        return np.nan
    t = rX * np.sqrt(df / (1 - rX**2))
    return float(2 * t_dist.sf(abs(t), df))
