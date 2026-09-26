import warnings

import numpy as np
import pytest

import pasta
from pasta.base import stats2nef, variogram_spec
from pasta.config import set_blocks
from pasta.distance_variable import load_distance_blockij_in_memory, prepare_distance_in_memory


@pytest.fixture
def maps():
    rng = np.random.default_rng(123)
    coord = rng.normal(size=(18, 3))
    x = rng.normal(size=18)
    y = 0.5 * x + rng.normal(size=18)
    return x, y, coord


def test_block_setting_and_variogram_spec():
    blocks = set_blocks(11, 4, 2)
    assert blocks.n_blocks == 3
    assert list(zip(blocks.block_i, blocks.block_j)) == [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    spec = variogram_spec(1, 10, 7, 0.7, 4, True)
    assert spec.h.shape == (7,)
    assert spec.L > 0 and spec.kernel_nh > 0


def test_packed_distance_uses_numpy_triu_order():
    N = 6
    D = np.zeros((N, N))
    upper = np.arange(1, N * (N - 1) // 2 + 1, dtype=float)
    D[np.triu_indices(N, 1)] = upper
    D += D.T
    blocks = set_blocks(N, 3)
    memory = prepare_distance_in_memory(upper, N)
    assert np.array_equal(load_distance_blockij_in_memory(memory, 0, 0, blocks), D[:3, :3])
    assert np.array_equal(load_distance_blockij_in_memory(memory, 0, 1, blocks), D[:3, 3:])


def test_distance_representations_and_parallel_are_equivalent(maps):
    x, y, coord = maps
    D = np.linalg.norm(coord[:, None, :] - coord[None, :, :], axis=-1)
    common = dict(M=7, block_size=5)
    euclidean = pasta.pasta(x, y, coord, **common)
    matrix = pasta.pasta(x, y, coord, D=D, **common)
    packed = pasta.pasta(x, y, coord, D=D[np.triu_indices(D.shape[0], 1)], **common)
    parallel = pasta.pasta(x, y, coord, n_workers=2, **common)
    for other in (matrix, packed, parallel):
        assert other.nef == pytest.approx(euclidean.nef, rel=1e-6, abs=1e-6)
        assert other.pef == pytest.approx(euclidean.pef, rel=1e-6, abs=1e-6)
    assert not hasattr(euclidean, "fc_para1")
    assert euclidean.model_x.parcel_b.shape == (1, 4)


def test_manual_parcels_and_invalid_masking(maps):
    x, y, coord = maps
    x[2] = np.nan
    labels = np.array([10] * 9 + [20] * 9)
    fit = pasta.pasta(x, y, coord, M=7, block_size=4, xparc=labels, yparc=labels)
    assert fit.n_parc.tolist() == [2, 2]
    assert fit.model_x.parc_idx.min() == 0
    assert fit.model_x.parc_idx.max() == 1
    assert fit.model_x.parc_idx.size == x.size - 1


def test_packed_distance_preserves_original_indices_after_masking(maps):
    x, y, coord = maps
    D = np.linalg.norm(coord[:, None, :] - coord[None, :, :], axis=-1)
    x[[1, 11]] = np.nan
    common = dict(M=7, block_size=4)
    matrix = pasta.pasta(x, y, coord, D=D, **common)
    packed = pasta.pasta(x, y, coord, D=D[np.triu_indices(D.shape[0], 1)], **common)
    assert packed.nef == pytest.approx(matrix.nef, rel=1e-6, abs=1e-6)


def test_auto_parcels_are_reproducible(maps):
    x, y, coord = maps
    options = dict(M=7, block_size=5, xparc="auto", yparc="auto", min_cluster_size=4, max_clusters=3, random_state=5)
    one = pasta.pasta(x, y, coord, **options)
    two = pasta.pasta(x, y, coord, **options)
    assert np.array_equal(one.model_x.parc_idx, two.model_x.parc_idx)
    assert np.array_equal(one.model_y.parc_idx, two.model_y.parc_idx)


def test_legacy_dense_api_is_deprecated(maps):
    x, y, coord = maps
    with pytest.warns(DeprecationWarning):
        result = pasta.effective_sample_size_estimation(x, y, coord, M=7)
    assert result[-2:] == (None, None)
    with pytest.warns(DeprecationWarning), pytest.raises(NotImplementedError):
        pasta.covariance_estimation(x, coord)


def test_stats2nef_invalid_covariance_product():
    nef, status = stats2nef(np.ones(3), np.ones(3), 1, 1, 1)
    assert not status
    assert not np.isfinite(nef)
