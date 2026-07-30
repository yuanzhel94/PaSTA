function pasta_fit = pasta(x, y, coord, varargin)
%pasta_fit = pasta(x, y, coord, varargin)
%This function run PaSTA in memory efficient way by evaluating distance
%and covariance matrices in blockwise manner

positive_int = @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x >= 1 && fix(x) == x;
p = inputParser;
addParameter(p, 'M', [], @(x) isempty(x) || positive_int(x));
addParameter(p, 'M_scale', 3, positive_int);
addParameter(p, 'qd', 0.7, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 1);
addParameter(p, 'nugget', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'kernel_scale', 4, @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0);
addParameter(p, 'block_size', 2000, @(x) isempty(x) || positive_int(x));
addParameter(p, 'n_workers', 1, positive_int);
addParameter(p, 'xparc', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)) || (ischar(x) && strcmp(x, 'auto')) || isstring(x) && isscalar(x) && x == "auto"); 
addParameter(p, 'yparc', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)) || (ischar(x) && strcmp(x, 'auto')) || isstring(x) && isscalar(x) && x == "auto");
addParameter(p, 'max_clusters', 10, @(x) isnumeric(x) && isscalar(x) && x > 0 && mod(x,1)==0);% only work when 'auto' is used
addParameter(p, 'min_clusters', 1, @(x) isnumeric(x) && isscalar(x) && x > 0 && mod(x,1)==0);% only work when 'auto' is used
addParameter(p, 'min_cluster_size', 500, @(x) isnumeric(x) && isscalar(x) && x > 0 && mod(x,1)==0); % only work when 'auto' is used
addParameter(p, 'D', [], @(x) isempty(x) || (isnumeric(x) && isreal(x) && ismatrix(x) && ~issparse(x)));
addParameter(p, 'distance_file', [], @(x) isempty(x) || is_v73_matfile(x));
addParameter(p, 'distance_variable', [], @(x) isempty(x) || ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'dim', [], @(x) isempty(x) || positive_int(x));
addParameter(p, 'random_state', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && fix(x) == x));
parse(p, varargin{:});

if ~isempty(p.Results.random_state)
    current_rng = rng;
    rng_cleanup = onCleanup(@() rng(current_rng));
    rng(p.Results.random_state);
end

distance_file = p.Results.distance_file;
distance_variable = p.Results.distance_variable;
D = p.Results.D;
use_distance_in_memory = ~isempty(D);
if isstring(distance_file)
    distance_file = char(distance_file);
end
if isstring(distance_variable)
    distance_variable = char(distance_variable);
end
if use_distance_in_memory
    load_distance = false;
else
    has_distance_file = ~isempty(distance_file);
    has_distance_variable = ~isempty(distance_variable);
    if xor(has_distance_file, has_distance_variable)
        error('distance_file and distance_variable must either both be empty or both be provided');
    end
    load_distance = has_distance_file && has_distance_variable;
end

xparc = p.Results.xparc;
yparc = p.Results.yparc;
x = x(:);
y = y(:);
N_in = length(x);
if length(y) ~= N_in
    error('x and y must be of same length');
end
if ~isnumeric(coord) || ~ismatrix(coord) || ~isreal(coord)
    error('coord must be a real numeric matrix');
elseif size(coord, 1) ~= N_in
    error('coord must have one row per observation');
end
valid = isfinite(x) & isfinite(y) & all(isfinite(coord), 2);

if isnumeric(xparc) && ~isempty(xparc)
    if numel(xparc) ~= N_in
        error('user-specified parcellation xparc must have one entry per original observation');
    end
    xparc = xparc(:);
    xparc = xparc(valid);
end
if isnumeric(yparc) && ~isempty(yparc)
    if numel(yparc) ~= N_in
        error('user-specified parcellation yparc must have one entry per original observation');
    end
    yparc = yparc(:);
    yparc = yparc(valid);
end
x = x(valid);
y = y(valid);
coord = coord(valid, :);
N = length(x);
[rX, p_naive] = corr(x, y);
std_x = std(x);
std_y = std(y);
if ~isfinite(std_x) || std_x <= 0
    error('x has zero or invalid variance, check input');
end
if ~isfinite(std_y) || std_y <= 0
    error('y has zero or invalid variance, check input.');
end
x = (x - mean(x)) ./ std_x;
y = (y - mean(y)) ./ std_y;
if isempty(p.Results.dim)
    dim = size(coord, 2);
else
    dim = double(p.Results.dim);
end
setting = pasta_setting(N, ...
    'M', p.Results.M, ...
    'M_scale', p.Results.M_scale, ...
    'qd', p.Results.qd, ...
    'nugget', p.Results.nugget, ...
    'kernel_scale', p.Results.kernel_scale, ...
    'block_size', p.Results.block_size, ...
    'n_workers', p.Results.n_workers, ...
    'xparc', xparc, ...
    'yparc', yparc, ...
    'max_clusters', p.Results.max_clusters, ...
    'min_clusters', p.Results.min_clusters, ...
    'min_cluster_size', p.Results.min_cluster_size);

if use_distance_in_memory
    distance_memory = prepare_distance_in_memory(D, N_in);
    distance_blocks = set_blocks(N_in, p.Results.block_size, p.Results.n_workers);
    distance_blocks.valid = valid(:);
    distance_blocks.in_idx_to_valid_idx = zeros(N_in, 1);
    distance_blocks.in_idx_to_valid_idx(valid) = (1:N)';
    [dmin, dmax] = distance_limits_with_blocks_in_memory(distance_memory, distance_blocks);
    spec = variogram_spec(dmin, dmax, setting.M, setting.qd, setting.kernel_scale, setting.nugget);
    [vx, vy] = global_variogram_with_blocks_in_memory(x, y, distance_memory, distance_blocks, spec);
    [global_bx, fx] = fit_stable_variogram(vx, spec.h, setting.nugget);
    [global_by, fy] = fit_stable_variogram(vy, spec.h, setting.nugget);
    global_bx = global_bx(:)';
    global_by = global_by(:)';
    n_points_from_range_x = NaN;
    n_points_from_range_y = NaN;
    use_range_x = ~isempty(setting.xparc);
    use_range_y = ~isempty(setting.yparc);
    ranges = zeros(use_range_x + use_range_y, 1);
    range_idx = 0;
    if use_range_x
        range_idx = range_idx + 1;
        ranges(range_idx) = global_bx(2) * 2.996 ^ (1 / global_bx(3));
    end
    if use_range_y
        range_idx = range_idx + 1;
        ranges(range_idx) = global_by(2) * 2.996 ^ (1 / global_by(3));
    end
    if ~isempty(ranges)
        counts = within_range_count_in_memory(ranges, distance_memory, distance_blocks);
        range_idx = 0;
        if use_range_x
            range_idx = range_idx + 1;
            n_points_from_range_x = counts(range_idx) / N - 1;
        end
        if use_range_y
            range_idx = range_idx + 1;
            n_points_from_range_y = counts(range_idx) / N - 1;
        end
    end
    [parc_idx_x, ~, n_parc_x] = get_parc(setting.xparc, coord, n_points_from_range_x, setting, 1);
    [parc_idx_y, ~, n_parc_y] = get_parc(setting.yparc, coord, n_points_from_range_y, setting, 2);
    [bs_x, bs_y] = fit_parcel_variogram_models_in_memory(x, y, distance_memory, parc_idx_x, n_parc_x, parc_idx_y, n_parc_y, global_bx, global_by, setting, distance_blocks);
    if n_parc_x > 1 || n_parc_y > 1
        [rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = prepare_nonstationary_stats_with_blocks_in_memory_distance(distance_memory, parc_idx_x, bs_x, parc_idx_y, bs_y, dim, distance_blocks);
    else
        [rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = prepare_stats_with_blocks_in_memory_distance(distance_memory, global_bx, global_by, distance_blocks);
    end

elseif load_distance
    distance_matfile = matfile(distance_file);
    variables_in_file = who(distance_matfile);
    if ~ismember(distance_variable, variables_in_file)
        error('Variable %s was not found in the distance MAT-file.', distance_variable);
    end
    distance_size = size(distance_matfile, distance_variable);
    if numel(distance_size) ~= 2 || distance_size(1) ~= distance_size(2)
        error('The precomputed distance variable must be square');
    end
    if distance_size(1) ~= N_in
        error('The precomputed distance matrix size must equal the number of inputs in x, including invalid values');
    end
    distance_blocks = set_blocks(N_in, p.Results.block_size, p.Results.n_workers);
    distance_blocks.valid = valid(:);
    distance_blocks.in_idx_to_valid_idx = zeros(N_in, 1);
    distance_blocks.in_idx_to_valid_idx(valid) = (1:N)';
    [dmin, dmax] = distance_limits_with_blocks_loaded(distance_matfile, distance_variable, distance_blocks);
    spec = variogram_spec(dmin, dmax, setting.M, setting.qd, setting.kernel_scale, setting.nugget);
    [vx, vy] = global_variogram_with_blocks_loaded(x, y, distance_matfile, distance_variable, distance_blocks, spec);
    [global_bx, fx] = fit_stable_variogram(vx, spec.h, setting.nugget);
    [global_by, fy] = fit_stable_variogram(vy, spec.h, setting.nugget);
    global_bx = global_bx(:)';
    global_by = global_by(:)';
    n_points_from_range_x = NaN;
    n_points_from_range_y = NaN;
    use_range_x = ~isempty(setting.xparc);
    use_range_y = ~isempty(setting.yparc);
    ranges = zeros(use_range_x + use_range_y, 1);
    range_idx = 0;
    if use_range_x
        range_idx = range_idx + 1;
        ranges(range_idx) = global_bx(2) * 2.996 ^ (1 / global_bx(3));
    end
    if use_range_y
        range_idx = range_idx + 1;
        ranges(range_idx) = global_by(2) * 2.996 ^ (1 / global_by(3));
    end
    if ~isempty(ranges)
        counts = within_range_count_loaded_distance(distance_matfile, distance_variable, ranges, distance_blocks);
        range_idx = 0;
        if use_range_x
            range_idx = range_idx + 1;
            n_points_from_range_x = counts(range_idx) / N - 1;
        end
        if use_range_y
            range_idx = range_idx + 1;
            n_points_from_range_y = counts(range_idx) / N - 1;
        end
    end
    [parc_idx_x, ~, n_parc_x] = get_parc(setting.xparc, coord, n_points_from_range_x, setting, 1);
    [parc_idx_y, ~, n_parc_y] = get_parc(setting.yparc, coord, n_points_from_range_y, setting, 2);
    [bs_x, bs_y] = fit_parcel_variogram_models_loaded(x, y, distance_matfile, distance_variable,...
        parc_idx_x, n_parc_x, parc_idx_y, n_parc_y, global_bx, global_by, setting, distance_blocks);
    if n_parc_x > 1 || n_parc_y > 1
        [rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = ...
            prepare_nonstationary_stats_with_blocks_loaded_distance(distance_matfile, ...
                distance_variable, ...
                parc_idx_x, ...
                bs_x, ...
                parc_idx_y, ...
                bs_y, ...
                dim, ...
                distance_blocks);
    else
        [rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = ...
            prepare_stats_with_blocks_loaded_distance(distance_matfile, ...
                distance_variable, ...
                global_bx, ...
                global_by, ...
                distance_blocks);
    end

% compute euclidean distance as you go   
else
    coord_norm2 = sum(coord .^ 2, 2);
    [dmin, dmax] = distance_limits_with_blocks_euclidean(coord, coord_norm2, setting);
    spec = variogram_spec(dmin, dmax, setting.M, setting.qd, setting.kernel_scale, setting.nugget);
    [vx, vy] = global_variogram_with_blocks_euclidean(x, y, coord, coord_norm2, setting, spec);
    [global_bx, fx] = fit_stable_variogram(vx, spec.h, setting.nugget);
    [global_by, fy] = fit_stable_variogram(vy, spec.h, setting.nugget);
    global_bx = global_bx(:)';
    global_by = global_by(:)';
    n_points_from_range_x = NaN;
    n_points_from_range_y = NaN;
    use_range_x = ~isempty(setting.xparc);
    use_range_y = ~isempty(setting.yparc);
    if use_range_x
        range_x = global_bx(2) * 2.996 ^ (1 / global_bx(3));
        count_x = within_range_count_euclidean(coord, coord_norm2, range_x, setting);
        n_points_from_range_x = count_x / N - 1;
    end
    if use_range_y
        range_y = global_by(2) * 2.996 ^ (1 / global_by(3));
        count_y = within_range_count_euclidean(coord, coord_norm2, range_y, setting);
        n_points_from_range_y = count_y / N - 1;
    end
    [parc_idx_x, ~, n_parc_x] = get_parc(setting.xparc, coord, n_points_from_range_x, setting, 1);
    [parc_idx_y, ~, n_parc_y] = get_parc(setting.yparc, coord, n_points_from_range_y, setting, 2);
    bs_x = variogram_fit_with_blocks_euclidean_all_parcels(x, coord, coord_norm2, parc_idx_x, n_parc_x, global_bx, setting);
    bs_y = variogram_fit_with_blocks_euclidean_all_parcels(y, coord, coord_norm2, parc_idx_y, n_parc_y, global_by, setting);
    if n_parc_x > 1 || n_parc_y > 1
        [rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = ...
            prepare_nonstationary_stats_with_blocks_euclidean(coord, ...
                coord_norm2, ...
                parc_idx_x, ...
                bs_x, ...
                parc_idx_y, ...
                bs_y, ...
                dim, ...
                setting);

    else
        [rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = ...
            prepare_stats_with_blocks_euclidean(coord, ...
                coord_norm2, ...
                global_bx, ...
                global_by, ...
                setting);
    end
end

[nef, run_status] = stats2nef(rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum);

if run_status
    pef = nef2p(rX, nef);
else
    pef = NaN;
end

pasta_fit = struct();
pasta_fit.pef = pef;
pasta_fit.rX = rX;
pasta_fit.nef = nef;
pasta_fit.run_status = run_status;
pasta_fit.n_parc = [n_parc_x, n_parc_y];
pasta_fit.p_naive = p_naive;

model_x = struct();
model_x.global_b = global_bx;
model_x.parcel_b = bs_x;
model_x.parc_idx = parc_idx_x;
model_x.n_parc = n_parc_x;
model_x.v = vx;
model_x.h = spec.h;
model_x.f = fx;

model_y = struct();
model_y.global_b = global_by;
model_y.parcel_b = bs_y;
model_y.parc_idx = parc_idx_y;
model_y.n_parc = n_parc_y;
model_y.v = vy;
model_y.h = spec.h;
model_y.f = fy;

pasta_fit.model_x = model_x;
pasta_fit.model_y = model_y;

end
