function [b, v, f] = variogram_fit_with_blocks_euclidean_single_parcel(x, coord, coord_norm2, global_b, setting)
%[b, v, f] = variogram_fit_with_blocks_euclidean_single_parcel(x, coord, coord_norm2, global_b, setting)
%This function estimate and fit the variogram for one parcel, by treating a parcel as a map,
%and using parcel_setting and parcel_spec for locally derived from the parcel data.

parcel_variance = var(x);
x = (x - mean(x)) ./ sqrt(parcel_variance);
N_parcel = length(x);
parcel_setting = pasta_setting(N_parcel, ...
    'block_size', min(setting.block_size, N_parcel), ...
    'n_workers', 1, ...
    'M_scale', setting.M_scale, ...
    'qd', setting.qd, ...
    'nugget', setting.nugget, ...
    'kernel_scale', setting.kernel_scale);

[dmin, dmax] = distance_limits_with_blocks_euclidean(coord, coord_norm2, parcel_setting);
parcel_spec = variogram_spec(dmin, dmax, parcel_setting.M, parcel_setting.qd, parcel_setting.kernel_scale, parcel_setting.nugget);

weight_sum = zeros(parcel_spec.M, 1);
weighted_diff2_sum = zeros(parcel_spec.M, 1);

for pair_idx = 1:parcel_setting.n_unique_block_pairs
    block_i = parcel_setting.block_i(pair_idx);
    block_j = parcel_setting.block_j(pair_idx);
    idx_i = parcel_setting.block_start(block_i):parcel_setting.block_end(block_i);
    idx_j = parcel_setting.block_start(block_j):parcel_setting.block_end(block_j);
    is_diagonal_block = parcel_setting.is_diagonal_block(pair_idx);
    d2ij = blockij_euclidean_d2(coord(idx_i,:), ...
        coord(idx_j,:), ...
        coord_norm2(idx_i), ...
        coord_norm2(idx_j), ...
        is_diagonal_block);
    dij = sqrt(d2ij);
    [weight_sum_ij, weighted_diff2_sum_ij] = blockij_variogram_accumulation_single_map(dij, ...
        x(idx_i), x(idx_j), parcel_spec, is_diagonal_block);
    weight_sum = weight_sum + weight_sum_ij;
    weighted_diff2_sum = weighted_diff2_sum + weighted_diff2_sum_ij;
end

v = nan(parcel_spec.M, 1);
valid_h = weight_sum > 0;
v(valid_h) = 0.5 .* weighted_diff2_sum(valid_h) ./ weight_sum(valid_h);

[b, f] = fit_stable_variogram_fixed_exponent(v, parcel_spec.h, global_b(3), parcel_setting.nugget);
b(1) = b(1) * parcel_variance;
b(4) = b(4) * parcel_variance;

end