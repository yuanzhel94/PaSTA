Getting Started
======================

Minimum data requirements for running PaSTA algorithms include the two spatial map observations :math:`x` and :math:`y`, and the spatial coordinates of observations :math:`coord`.
When :math:`coord` is not available, the algorithm can also run with the pairwise distance matrix between observations :math:`D` (see **Advanced examples** for more details).

:math:`x`: ndarray (N,) in python; vector (N,1) in MATLAB. Can contain NaN and Inf.

:math:`y`: ndarray (N,) in python; vector (N,1) in MATLAB. Can contain NaN and Inf.

:math:`coord`: ndarray (N,dim) in python; vector (N,dim) in MATLAB, where dim is 3 for 3D data.

In this page, we provide guidance to run **PaSTA** and **PaSTA-NS** using the default settings, in MATLAB and Python. 
This default setting should be sufficient if you aim to:

1. Run PaSTA, or
2. Run PaSTA-NS for large covareage (e.g., whole-brain) maps on a resolution comparable to fsavearge5 10k vertices map.

Otherwise, I recommend see  **Advanced examples** for more details on appropriate use. These advanced examples include but not limited to:

1. Run PaSTA-NS for small data patches (see **'Controlling for parcels in PaSTA-NS'** in **Advanced examples**).
2. Run PaSTA-NS using user-defined parcels (see **'Controlling for parcels in PaSTA-NS'** in **Advanced examples**).
3. When :math:`coord` is missing (see **'Run with distance matrix'** in **Advanced examples**).
4. Run pairwise associations between a large number of spatial maps (see **'Large-scale pairwise evaluation'** in **Advanced examples**).
5. In the presence of spatial trends (see **'Impact of spatial trends'** in **Advanced examples**).

PaSTA
----------------------------
.. tabs:: lang

    .. code-tab:: MATLAB

        % take 70-80s to run (Apple Silicon M1 Pro) on fsaverage5 10k cortical map
        % pef - significance p-value
        % rX - Pearson correlation coefficient
        % nef - effective sample size
        [pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2] = effective_sample_size_estimation(x,y,coord);

    .. code-tab:: python

        # take 3.5min to run (Apple Silicon M1 Pro) on fsaverage5 10k cortical map
        # pef - significance p-value
        # rX - Pearson correlation coefficient
        # nef - effective sample size
        import pasta
        pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2 = pasta.effective_sample_size_estimation(x, y, coord)


PaSTA-NS
----------------------------
.. tabs:: lang

    .. code-tab:: MATLAB

        % PaSTA-NS with data-driven parcellation
        % take 70-80s to run (Apple Silicon M1 Pro) on fsaverage5 10k cortical map
        % pef - significance p-value
        % rX - Pearson correlation coefficient
        % nef - effective sample size
        [pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2] = effective_sample_size_estimation(x,y,coord,'xparc','auto','yparc','auto');

    .. code-tab:: python

        # PaSTA-NS with data-driven parcellation
        # take 3.5min to run (Apple Silicon M1 Pro) on fsaverage5 10k cortical map
        # pef - significance p-value
        # rX - Pearson correlation coefficient
        # nef - effective sample size
        import pasta
        pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2 = pasta.effective_sample_size_estimation(x, y, coord, xparc='auto', yparc='auto')
