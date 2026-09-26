"""Deprecated compatibility API for the pre-1.0 dense Python package."""
from __future__ import annotations

import warnings
from .pasta import pasta


def effective_sample_size_estimation(x, y, coord=None, D=None, dim=None, M=None, qd=0.7, xparc=None, yparc=None, max_clusters=10, min_cluster_size=500, min_clusters=1, M_cluster=None, nugget=True):
    """Deprecated wrapper; covariance outputs are always ``None``."""
    warnings.warn("effective_sample_size_estimation is deprecated; use pasta.pasta, which returns PaSTAFit and never materializes covariance matrices", DeprecationWarning, stacklevel=2)
    if coord is None:
        raise ValueError("coord is required by the MATLAB-aligned pasta API")
    if M_cluster is not None:
        warnings.warn("M_cluster is ignored; MATLAB PaSTA uses M_scale for parcel variograms", DeprecationWarning, stacklevel=2)
    fit = pasta(x, y, coord, M=M, qd=qd, xparc=xparc, yparc=yparc, max_clusters=max_clusters, min_clusters=min_clusters, min_cluster_size=min_cluster_size, D=D, dim=dim, nugget=nugget)
    return fit.pef, fit.rX, fit.nef, fit.run_status, fit.n_parc, fit.p_naive, None, None


def _dense_removed(name):
    warnings.warn(f"{name} is deprecated because it requires dense covariance matrices; PaSTA 1.0 never materializes N-by-N covariance matrices", DeprecationWarning, stacklevel=2)
    raise NotImplementedError(f"{name} was removed; use pasta.pasta for blockwise inference")


def covariance_estimation(*args, **kwargs): _dense_removed("covariance_estimation")
def cov2nef(*args, **kwargs): _dense_removed("cov2nef")
def estimate_variogram(*args, **kwargs): _dense_removed("estimate_variogram")
def fit_variogram(*args, **kwargs): _dense_removed("fit_variogram")
def fit_variogram_fixed_exponent(*args, **kwargs): _dense_removed("fit_variogram_fixed_exponent")
def fit_covariance_blocks(*args, **kwargs): _dense_removed("fit_covariance_blocks")
def process_convolution_crossblocks(*args, **kwargs): _dense_removed("process_convolution_crossblocks")
