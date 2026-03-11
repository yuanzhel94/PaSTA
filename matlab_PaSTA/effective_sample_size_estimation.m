function [pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2] = effective_sample_size_estimation(x,y,coord,varargin)
% This function computes the significance for spatial association using PaSTA
% Inputs:
%   x: frist spatial data vector of length N, can contain NaN and Inf
%   y: second spatial data vector of length N, can contain NaN and Inf
%   coord: spatial coordinates of the data in x and y of shape (N, dim)
% Optional Inputs:
%   M: number of lag distances to evaluate when estimating the global
%       variogram, default to 3*sqrt(N) if not provided
%   qd: factor determines the largest distance to evaluate in variogram,
    %   float number (0,1]. Default to 0.7. Largest distance evaluated is
    %   computed as qd*max(D(:)), i.e., the product of qd and the maximum
    %   distance between data points.
%   xparc: either [], 'auto', or int vector of length N where each value represents a parcellation.
    %   Determines the way to parcellate data x and compute nonstationary covariance. Default to []
    %   that assumes stationary autocorrelation and does not parcellate (i.e.,
    %   PaSTA). if 'auto', run PaSTA-NS by determining the parcel using a
    %   data-driven spatial clustering. if vector of length N, parcellate data
    %   using user-specified parcellation to estimate nonstationarity.
%   yparc: same as xparc but determines the way to parcellate data y
%   max_clusters: the maximum number of parcels allowed if data-driven
%       parcellation. Default to 10. Only work when xparc or yparc == 'auto'.
%   min_clusters: the minimum number of parcels allowed if data-driven
    %   parcellation. Default to 1. This is used in our paper to investigate
    %   the impact of mandatory parcellation when data-driven approach
    %   suggest autocorrelation is too strong to allow data parcellation. In
    %   practice, always keep to 1 that is the default. Only work when xparc or yparc == 'auto'.
%   min_cluster_size: the average number of samples in each parcel if
%       data-driven parcellation. Only work when xparc or yparc == 'auto'.
%   crossblock: ways to infer between-parcel covariance. Although 'conv'
    %   (process convolution) and 'independent' (assuming spatial independence
    %   between parcels) are implemented, should always keep to 'conv', which
    %   is the default. Only work in PaSTA-NS, i.e., number of parcels are > 1.
%   nugget: bool indicator of whether use nugget in variogram models.
    %   Should always keep to true (default) for better variogram fitting.
%   M_scale: the scaling factor to determine the number of lag distances
    %   evaluated in variograms if M not specified. default to 3 which provides
    %   a good FPR control in both short and long range autocorrelation
    %   settings.
% Outputs:
%   pef: statistical significance p-value derived using PaSTA/PaSTA-NS
%   rX: pearson correlation coefficient
%   nef: effective sample size estiamted with PaSTA/PaSTA-NS
%   run_status: bool indicator whether PaSTA/PaSTA-NS estimation is valid, 1
    %   indicates success, 0 indicates unsuccess which can happen when data is 
    %   too smooth and nef < 2
%   n_parc: [nx_parc, n_yparc] that indicates the number of parcels for
    %   each map is used in PaSTA-NS
%   p_naive: statistical significance not accounting for spatial
    %   autocorrelation and assuming spatial independence
%   fc_para1: estimated covariance matrix for spatial map x
%   fc_para2: estimated covariance matrix for spatial map y

p = inputParser;
addParameter(p, 'M', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v > 0 && mod(v,1)==0));
addParameter(p, 'qd', 0.7, @(v) isnumeric(v) && v > 0 && v <= 1);
addParameter(p, 'xparc', [], @(v) (isnumeric(v) && isvector(v)) || (ischar(v) && strcmp(v, 'auto')));
addParameter(p, 'yparc', [], @(v) (isnumeric(v) && isvector(v)) || (ischar(v) && strcmp(v, 'auto')));
addParameter(p, 'max_clusters', 10, @(v) isnumeric(v) && isscalar(v) && v > 0 && mod(v,1)==0);% only work when 'auto' is used
addParameter(p, 'min_clusters', 1, @(v) isnumeric(v) && isscalar(v) && v > 0 && mod(v,1)==0);% only work when 'auto' is used
addParameter(p, 'min_cluster_size', 500, @(v) isnumeric(v) && isscalar(v) && v > 0 && mod(v,1)==0); % only work when 'auto' is used
addParameter(p, 'crossblock', 'conv', @(v) ismember(v, {'independent', 'conv'})); % conv/indepdent using convolution process/indepdence for cross-parcs/cross-blocks
addParameter(p, 'nugget', true, @(v) islogical(v)); % use nugget
addParameter(p, 'M_scale', 3, @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v > 0 && mod(v,1)==0));
parse(p, varargin{:});

fittype = 'stable'; % fit varigorams with stable model

M = p.Results.M;
M_scale = p.Results.M_scale;
if isempty(M)
    M = M_scale * ceil(sqrt(length(x))); % if number of bins (M) not specified, using the default as M_scale * sqrt(len(data))
end
qd = p.Results.qd;
xparc = p.Results.xparc;
yparc = p.Results.yparc;
crossblock = p.Results.crossblock;
nugget = p.Results.nugget;
max_clusters = p.Results.max_clusters;
min_clusters = p.Results.min_clusters;
min_cluster_size = p.Results.min_cluster_size;

% remove nan and inf data points before computation
valid = (isfinite(x) & isfinite(y));
x = x(valid); y = y(valid); 
coord = coord(valid,:);
D = squareform(pdist(coord));
dim = size(coord,2);

if (~isempty(xparc)) & (~strcmp(xparc, 'auto'))
    xparc = xparc(valid);
end
if (~isempty(yparc)) & (~strcmp(yparc, 'auto'))
    yparc = yparc(valid);
end

x = zscore(x); y = zscore(y);
PrecomputedVariance = [];

% now run the stationary stable model without parcellation, i.e., PaSTA
[v1,h1] = estimate_variogram(D, x, M, qd);
[v2,h2] = estimate_variogram(D, y, M, qd);
[c_para1, b1]=fit_variogram(h1',v1,D,fittype,PrecomputedVariance,nugget);
[c_para2, b2]=fit_variogram(h2',v2,D,fittype,PrecomputedVariance,nugget);

% parcelating x and y if needed, fc_para is returned as c_para if parcel is
% not needed, otherwise zeros matrix
[xparc, xn_parc, xunique_parcs, fc_para1] = parc_data(xparc, c_para1, b1, D, coord, max_clusters, min_clusters, min_cluster_size, 1);
[yparc, yn_parc, yunique_parcs, fc_para2] = parc_data(yparc, c_para2, b2, D, coord, max_clusters, min_clusters, min_cluster_size, 2);

% estimate the covariance matrix by parcellation if n_parc > 1, i.e., PaSTA-NS
if xn_parc > 1
    [fc_para1, bs1] = parcellated_cov(x, xparc, xunique_parcs, D, dim, b1, fc_para1, qd, M_scale, nugget, crossblock);
end
if yn_parc > 1
    [fc_para2, bs2] = parcellated_cov(y, yparc, yunique_parcs, D, dim, b2, fc_para2, qd, M_scale, nugget, crossblock);
end

[nef, run_status] = covs2nef(fc_para1, fc_para2);
n_parc = [xn_parc, yn_parc];
[rX, p_naive] = corr(x,y);
if run_status
    pef = nef2p(rX, nef);
else
    pef = NaN;
end

end
