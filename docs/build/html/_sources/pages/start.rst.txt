Getting Started
===============

Brain-PaSTA follows the MATLAB ``pasta`` entry point. Supply two spatial maps
and an ``(N, dim)`` coordinate array. Invalid rows are removed consistently
from both maps and the coordinates.

PaSTA
-----

.. code-block:: python

    import pasta

    fit = pasta.pasta(x, y, coord)
    print(fit.pef, fit.rX, fit.nef)

``fit`` is a ``PaSTAFit`` object with ``pef``, ``rX``, ``nef``,
``run_status``, ``n_parc``, ``p_naive``, ``model_x``, and ``model_y`` fields.
No covariance matrix is constructed or returned.

PaSTA-NS
--------

.. code-block:: python

    fit = pasta.pasta(x, y, coord, xparc="auto", yparc="auto", random_state=0)

Automatic parcels are generated from coordinates with K-means. User-defined
parcel labels are also accepted through ``xparc`` and ``yparc``.

Precomputed distances
---------------------

Provide ``D`` while still supplying coordinates. ``D`` can be a full ``(N, N)``
distance matrix or a NumPy-order strict-upper-triangle vector:

.. code-block:: python

    fit = pasta.pasta(x, y, coord, D=D)
    fit = pasta.pasta(x, y, coord, D=D[np.triu_indices(D.shape[0], k=1)])

The packed-vector ordering is NumPy's row-major ``triu_indices`` ordering, not
MATLAB's column-major linear-index ordering.
