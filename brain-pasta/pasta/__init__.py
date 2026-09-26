"""Brain-PaSTA's MATLAB-aligned Python interface."""

from .pasta import pasta
from .types import PaSTAFit, PaSTAModel
from .base import nef2p
from .legacy import covariance_estimation, cov2nef, effective_sample_size_estimation

__all__ = ["pasta", "PaSTAFit", "PaSTAModel", "nef2p", "effective_sample_size_estimation", "covariance_estimation", "cov2nef"]
