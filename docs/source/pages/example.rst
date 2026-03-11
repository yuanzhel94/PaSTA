Advanced examples
==================

This page demonstrates advanced examples of using PaSTA.

Controlling for parcels in PaSTA-NS
------------------------------------------------

Several arguments in **pasta.effective_sample_size_estimation()** and **pasta.covariance_estimation()** can be used to control the behavior of parcellations used in PaSTA-NS, including:
*x_parc* (*y_parc*), *max_clusters*, *min_clusters*, and *min_cluster_size*.

*x_parc* (and *y_parc*) determines whether PaSTA or PaSTA-NS is used, and also how data *x* (and *y*) are pacellated in PaSTA-NS. They can be one of the following 3 scenarios:

1. :math:`None`: PaSTA will be used.
2. String :math:`'auto'`: PaSTA-NS will be used, with parcellation in map is determined by data-driven spatial clustering, conditioned on parameters including: *max_clusters*, *min_clusters*, and *min_cluster_size*
3. int ndarray :math:`(N,)` with :math:`k` unique integers ranging from :math:`0` to :math:`k-1`, each indicate a parcellation: PaSTA-NS will be used, and parcellation is determined as use-specified. A warning will be raised if the number of user-specified parcellation is larger than the number of data-driven determined number conditioned on parameters including: *max_clusters*, *min_clusters*, and *min_cluster_size*.

*max_clusters (default 10)* and *min_cluster_size (default 500)* jointly determines the largest number of parcels can be used in data-driven parcellation with spatial clustering. The maximum number is the smaller of *max_clusters* or :math:`floor(N / {min\_cluster\_size})`. The *min_cluster_size* is included such that the average size (i.e., number of observations) of parcels is not too small to get accurate and reliable variogram estimation.

*min_clusters (default 1)* determines the smallest number of parcels can be used, which should be always set to 1 in practice. This argument is used in our study to investigate the impact of forced parcellation.

.. note::

    Use *x_parc* and *y_parc* to manually set the choice of parcellations in PaSTA-NS

.. warning::

    May need to change the input of *max_clusters* and *min_cluster_size* when apply PaSTA-NS with data-driven spatial clustering for small data patches, depending on your data.


Run with distance matrix (python only)
------------------------------------------------

PaSTA can be run when the spatial coordinates of observations are missing but the pairwise distances :math:`D` is known.

PaSTA-NS with data-driven spatial clustering can be run when the spatial coordinates of observations are missing but the pairwise distances :math:`D` and the spatial dimension :math:`dim` is known, where clusters are determined using K-Medoids clustering implemented with Scikit-learn-extra.

Examples below:

.. code-tab:: python

    import pasta
    # PaSTA
    pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2 = pasta.effective_sample_size_estimation(x, y, D=D)
    # PaSTA-NS
    pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2 = pasta.effective_sample_size_estimation(x, y, D=D, dim=dim, xparc='auto', yparc='auto')


Large-scale pairwise evaluation
-----------------------------------------------

The **pasta.effective_sample_size_estimation()** funciton internally estiamtes the covariance matrices of spatial maps, which can increase computational cost when evaluating pairwise associations between a large number of spatial maps because of repeatitive covariance estimation.
**pasta.covariance_estimation()** can be used in this case to improve computational efficiency, by decoupling covariance estimation and significance inference. In brief, the covariance matrices estimated can be saved and reloaded for later inference of statistical significance.

Below is an example of running this decoupled pipeline of PaSTA, to infer significance of associations between 3 different spatial maps.

Step 1: covariance estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. tabs:: lang

    .. code-tab:: MATLAB

        spatial_maps = {x, y, z};
        mapnames = {'x', 'y', 'z'};
        cov_path = '.' % save covariance matrices in the current path
        for i = 1:length(spatial_maps)
            mapi = spatial_maps{i};
            % computed covariance matrix may contain NaN if the spatial data contains NaN or Inf
            covmat = covariance_estimation(x,coord) % change function parameters for PaSTA-NS and other settings
            save(fullfile(cov_path, sprintf('%s_cov.mat', mapnames{i})), 'covmat')

    .. code-tab:: python

        import pasta
        import numpy as np
        import os

        spatial_maps = [x, y, z]
        mapnames = ['x', 'y', 'z']
        cov_path = '.' # save covariance matrix in the current dir
        for _, (mapi, namei) in enumerate(zip(spatial_maps, mapnames)):
            covmat = pasta.covariance_estimation(mapi, coord) # change function parameters for PaSTA-NS and other settings
            np.save(os.path.join(cov_path, f'{namei}_cov.npy'), covmat)
        
Step 2: infer significance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. tabs:: lang

    .. code-tab:: MATLAB

        spatial_maps = {x, y, z};
        mapnames = {'x', 'y', 'z'};
        cov_path = '.'
        ps = zeros(length(spatial_maps)) % initiate for saving
        rXs = zeros(length(spatial_maps)) % initiate for saving
        for i = 1:length(spatial_maps)-1
            fpr j = i+1:length(spatial_maps)
                mapi = spatial_maps{i};
                mapj = spatial_maps{j};
                covi = load(fullfile(cov_path, sprintf('%s_cov.mat', mapnames{i})));
                covj = load(fullfile(cov_path, sprintf('%s_cov.mat', mapnames{j})));
                % because data mapi and mapj, and covariance matrices covi and covj may contain NaN and Inf, remove these data points
                valid = isfinite(mapi) & isfinite(mapj);
                covi = covi(valid, valid);
                covj = covj(valid, valid);
                [rXs(i, j), ~] = corr(mapi(valid), mapj(valid));
                [nef, run_status] = covs2nef(covi,covj)
                ps(i,j) = nef2p(rXs(i,j), nef)

    .. code-tab:: python

        import pasta
        import numpy as np
        import os
        from scipy.stats import pearsonr

        spatial_maps = [x, y, z]
        mapnames = ['x', 'y', 'z']
        savepath = '.'
        n_maps = len(spatial_maps)
        rXs = np.zeros((n_maps, n_maps)) # initiate for saving
        ps = np.zeros((n_maps, n_maps)) # initiate for saving
        for i in np.arange(n_maps-1):
            for j in np.arange(i+1, n_maps):
                mapi = spatial_maps[i]
                namei = mapnames[i]
                mapj = spatial_maps[j]
                namej = spatial_maps[j]
                covi = np.load(os.path.join(cov_path, f'{namei}_cov.npy'))
                covj = np.load(os.path.join(cov_path, f'{namej}_cov.npy'))
                valid = np.logical_and(np.isfinite(mapi), np.isfinite(mapj))
                rX, p_naive = pearsonr(x[valid], y[valid])
                nef = pasta.cov2nef(covi[np.ix_(valid,valid)],covj[np.ix_(valid,valid)])
                pef = pasta.nef2p(rX, nef) if nef>2 else np.nan 
                rXs[i,j], ps[i,j] = rX, pef

Impact of spatial trends
-----------------------------------------------
Coming soon.
