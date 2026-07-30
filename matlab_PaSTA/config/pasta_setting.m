function setting = pasta_setting(N, varargin)
% setting = pasta_setting(N, varargin)
% This function specifies the computational setting for PaSTA
%
% Inputs:
%   N: Number of valid spatial observations.
%
% Optional Inputs:
%   block_size: size of block chunks in computation for memory saving, default N.
%   n_workers: number of parallel workers ,default 1.
%   M: number of lag distances to evaluate when estimating the global
%       variogram, default to M_scale*sqrt(N) if not provided
%   M_scale: scaling factor for M when M is determined from N, default to 3
%   qd: factor determines the largest distance to evaluate in variogram,
%      float number (0,1]. Default to 0.7. Largest distance evaluated is
%      computed as qd*max(D(:)), i.e., the product of qd and the maximum
%      distance between data points.
%   nugget: bool indicator of using nugget in variogram model, default true.
%       Should always keep to true (default) for better variogram fitting.
%   xparc: either [], 'auto', or int vector of length N where each value represents a parcellation.
%       Determines the way to parcellate data x and compute nonstationary
%       covariance. Default to [] that assumes stationary autocorrelation and does not parcellate (i.e.,
%       PaSTA). if 'auto', run PaSTA-NS by determining the parcel using a
%       data-driven spatial clustering. if vector of length N, parcellate data
%       using user-specified parcellation to estimate nonstationarity.
%   max_clusters: the maximum number of parcels allowed if data-driven
%       parcellation. Default to 10. Only work when xparc or yparc == 'auto'.
%   min_clusters: the minimum number of parcels allowed if data-driven
%       parcellation. Default to 1. This is used in our paper to investigate
%       the impact of mandatory parcellation when data-driven approach
%       suggest autocorrelation is too strong to allow data parcellation. In
%       practice, always keep to 1 that is the default. Only work when xparc or yparc == 'auto'.
%   min_cluster_size: the average number of samples in each parcel if
%       data-driven parcellation. Only work when xparc or yparc == 'auto'
%   kernel_scale: the scaling factor to determine the truncation window for
%       variogram estimation, efault to 4, i.e., 4 standard deviation.
%
%
% Outputs:
%   setting: setting struct specifying PaSTA settings, including
%       N, block_size, n_workers, parallel, stationarity, n_blocks, 
%       block_start, block_end, n_unique_block_pairs, block_i, block_j,
%       is_diagonal_block
%

p = inputParser;
addRequired(p, 'N');
addParameter(p, 'block_size', []);
addParameter(p, 'n_workers', 1);
addParameter(p, 'M', []);
addParameter(p, 'M_scale', 3);
addParameter(p, 'qd', 0.7);
addParameter(p, 'nugget', true);
addParameter(p, 'xparc', []); 
addParameter(p, 'yparc', []); 
addParameter(p, 'max_clusters', 10);% only work when 'auto' is used
addParameter(p, 'min_clusters', 1);% only work when 'auto' is used
addParameter(p, 'min_cluster_size', 500); % only work when 'auto' is used
addParameter(p, 'kernel_scale', 4);
addParameter(p, 'distance_file', []);
addParameter(p, 'distance_variable', []);
addParameter(p, 'valid_distance_entry_idx', []);
parse(p, N, varargin{:});

setting=struct();
setting.N = double(p.Results.N);
setting.M_scale = double(p.Results.M_scale);
setting.qd = double(p.Results.qd);
setting.kernel_scale = double(p.Results.kernel_scale);
setting.nugget = p.Results.nugget;

setting.xparc = p.Results.xparc;
setting.yparc = p.Results.yparc;

setting.max_clusters = double(p.Results.max_clusters);
setting.min_clusters = double(p.Results.min_clusters);
setting.min_cluster_size = double(p.Results.min_cluster_size);

if isempty(p.Results.M)
    setting.M = setting.M_scale * ceil(sqrt(setting.N));
else
    setting.M = double(p.Results.M);
end

if isempty(p.Results.block_size)
    setting.block_size = setting.N;
else
    setting.block_size = min(double(p.Results.block_size), setting.N);
end

setting.n_workers = double(p.Results.n_workers);

setting.distance_file = p.Results.distance_file;
setting.distance_variable = p.Results.distance_variable;
setting.load_distance = ~isempty(setting.distance_file) && ~isempty(setting.distance_variable);
setting.valid_distance_entry_idx = p.Results.valid_distance_entry_idx;

setting.n_blocks = ceil(setting.N / setting.block_size);
setting.block_start = (0:(setting.n_blocks - 1))' * setting.block_size + 1;
setting.block_end = min(setting.block_start + setting.block_size - 1, setting.N);
setting.block_length = setting.block_end - setting.block_start + 1;

setting.n_unique_block_pairs = setting.n_blocks * (setting.n_blocks + 1) / 2;
setting.block_i = nan(setting.n_unique_block_pairs, 1);
setting.block_j = nan(setting.n_unique_block_pairs, 1);
setting.parallel = (setting.n_workers > 1) && (setting.n_unique_block_pairs > 1);

pointer_idx = 1;
for i = 1:setting.n_blocks
    nj = setting.n_blocks - i + 1; % all j >= i
    idx = pointer_idx : (pointer_idx + nj -1);
    setting.block_i(idx) = i;
    setting.block_j(idx) = i:setting.n_blocks;
    pointer_idx = pointer_idx + nj;
end
setting.is_diagonal_block = setting.block_i == setting.block_j;

end

