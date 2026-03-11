.. PaSTA documentation master file, created by
   sphinx-quickstart on Sun Feb  1 20:49:02 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PaSTA documentation
=====================

PaSTA (Parametric Spatial Test for Associations) presents an analytical approach that infers significance of spatial associaiton tests between spatial maps by estimating the effective degrees of freedom of the data. 
While the development and validation focused neuroimaging data, the approach also works in the context of broader spatial data.

Algorithms
---------------------
PaSTA (assuming spatial stationarity/homogeneity)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Estimate the empirical variogram and fit a stable variogram model for each autocorrelated map :math:`x` and :math:`y`.
2. Compute the covariance matrices :math:`C_{x}` and :math:`C_{y}` based on the variogram.
3. Compute the effective sample size :math:`N_{\mathrm{ef}}` and effective degrees of freedom (:math:`N_{\mathrm{ef}} - 2`) for association test based on Dutilleul's derivation (Dutilleul, 1993):

.. math::

   N_{\mathrm{ef}} =
   \frac{\operatorname{tr}(B C_x B C_y)}
        {\operatorname{tr}(B C_x)\operatorname{tr}(B C_y)}
   + 1

4. Compute the statistical significance p-value by referencing the test statistic to its theoretical distribution based on the effective degrees of freedom


PaSTA-NS (can account for nonstationarity)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Estimate the empirical variogram and fit a stable variogram model for each autocorrelated map :math:`x` and :math:`y`.
2. Determine number of parcels based on the global strengths of autocorrelation (i.e., based on the :math:`range` parameter in the global variogram model).
3. Parcellate data using spatial clustering.
4. Estimate and fit a separate variogram model for each parcel.
5. Compute the nonstationary covariance matrice :math:`C_{x}` and :math:`C_{y}` from parcels based on process convolution (Paciorek & Schervish, 2006):

.. math::

   C^{NS}(i,j)
   =
   |\Sigma_i|^{1/4}
   |\Sigma_j|^{1/4}
   \left|\frac{\Sigma_i+\Sigma_j}{2}\right|^{-1/2}
   \;
   C^{S}\!\left(\sqrt{Q_{ij}}\right)

.. math::

   Q_{ij}
   =
   (s_i - s_j)^T
   \left(\frac{\Sigma_i+\Sigma_j}{2}\right)^{-1}
   (s_i - s_j)

where :math:`\Sigma_i` and :math:`\Sigma_j` are the covariance matrices of the
Gaussian kernel centered at spatial locations :math:`i` and :math:`j`,
respectively. :math:`C^{S}` denotes the stationary covariance function.

6. Compute the effective sample size and effective degrees of freedom for association test based on Dutillel's derivation using nonstationary covariance matrices.
7. Compute the statistical significance p-value by referencing the test statistic to its theoretical distribution based on the effective degrees of freedom


.. toctree::
   :maxdepth: 4
   :caption: Contents:
   
   pages/install
   pages/start
   pages/example
   api/docwrap
   pages/references
