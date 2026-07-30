function covmat = covariance_estimation(x,coord,varargin)
%This function estimates the covariance matrix for a given spatial map.
%This is particularly useful when estimating pairwise association between a
%large number of spatial maps, where the function
%effectie_sample_size_estimation can be redundant and repeatatively
%computes the covariance matrices for the same map. By saving the computed
%covariane matrix and load for inference involving the map, it reduces the
%computational cost in this setting.
% Inputs:
%   mostly similar to effective_sample_size_estimation but now does not
    %   need to provide 'y' and 'yparc' because covariance estimation only
    %   involve one map but significance inference require map pairs
% Outputs:
%   covmat: estimated covariance matrix of shape (N,N) for the input vector x of
    %   length N. Note this may contain NaN values for rows and columns
    %   corresponding to invalid data in x. When loading this covmat for
    %   significance inference, please remove the rows and columns
    %   corresponding to NaN values in the map pairs of interest.

p = inputParser;
addParameter(p, 'M', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v > 0 && mod(v,1)==0));
addParameter(p, 'qd', 0.7, @(v) isnumeric(v) && v > 0 && v <= 1);
addParameter(p, 'xparc', [], @(v) (isnumeric(v) && isvector(v)) || (ischar(v) && strcmp(v, 'auto')));
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
crossblock = p.Results.crossblock;
nugget = p.Results.nugget;
max_clusters = p.Results.max_clusters;
min_clusters = p.Results.min_clusters;
min_cluster_size = p.Results.min_cluster_size;

valid = isfinite(x);
nx = length(x);
covmat = NaN(nx);
x = x(valid);
coord = coord(valid,:);
D = squareform(pdist(coord));
dim = size(coord,2);

if (~isempty(xparc)) & (~strcmp(xparc, 'auto'))
    xparc = xparc(valid);
end
x = zscore(x);
PrecomputedVariance = [];

% now run the stationary stable model without parcellation
[v1,h1] = estimate_variogram(D, x, M, qd);
[c_para1, b1]=fit_variogram(h1',v1,D,fittype,PrecomputedVariance,nugget);

[xparc, xn_parc, xunique_parcs, fc_para1] = parc_data(xparc, c_para1, b1, D, coord, max_clusters, min_clusters, min_cluster_size, 1);

% estimate the covariance matrix by parcellation if n_parc > 1
if xn_parc > 1
    [fc_para1, bs1] = parcellated_cov(x, xparc, xunique_parcs, D, dim, b1, fc_para1, qd, M_scale, nugget, crossblock);
end
covmat(valid, valid) = fc_para1;
end