Advanced examples
=================

Controlling PaSTA-NS parcels
----------------------------

``xparc`` and ``yparc`` independently select stationary PaSTA (``None``),
data-driven PaSTA-NS (``"auto"``), or user-defined parcellation labels.

.. code-block:: python

    fit = pasta.pasta(
        x, y, coord,
        xparc="auto", yparc="auto",
        max_clusters=10,
        min_clusters=1,
        min_cluster_size=500,
        random_state=0,
    )

For a known parcellation, labels may be any one-dimensional values; PaSTA
normalizes them to zero-based indices in ``fit.model_x.parc_idx`` and
``fit.model_y.parc_idx``.

Memory and parallelism
----------------------

PaSTA allocates only an individual distance/covariance block at a time. Tune
the temporary-memory/performance trade-off with ``block_size``. Set
``n_workers`` above one to evaluate independent block pairs in parallel;
reductions remain deterministic.

.. code-block:: python

    fit = pasta.pasta(x, y, coord, block_size=2000, n_workers=4)

The package deliberately does not offer dense covariance estimation or a
distance-only clustering mode. Legacy dense functions emit deprecation
warnings and cannot allocate covariance matrices.
