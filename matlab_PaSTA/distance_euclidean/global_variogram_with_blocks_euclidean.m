function [vx, vy] = global_variogram_with_blocks_euclidean(x, y, coord, coord_norm2, setting, spec)
%[vx, vy] = global_variogram_with_blocks_euclidean(x, y, coord, coord_norm, setting, variogram_spec)
%This funciton computes the global variogram by aggregating blocks based on
%pasta setting and variogram specifications

if setting.parallel
    weight_sum_per_block = zeros(spec.M, setting.n_unique_block_pairs);
    weighted_diff2_sum_x_per_block = zeros(spec.M, setting.n_unique_block_pairs);
    weighted_diff2_sum_y_per_block = zeros(spec.M, setting.n_unique_block_pairs);
    
    parfor (pair_idx = 1:setting.n_unique_block_pairs, setting.n_workers)
        block_i = setting.block_i(pair_idx);
        block_j = setting.block_j(pair_idx);
        idx_i = setting.block_start(block_i):setting.block_end(block_i);
        idx_j = setting.block_start(block_j):setting.block_end(block_j);
        is_diagonal_block = setting.is_diagonal_block(pair_idx);
        d2 = blockij_euclidean_d2(coord(idx_i,:), ...
            coord(idx_j,:), ...
            coord_norm2(idx_i), ...
            coord_norm2(idx_j), ...
            is_diagonal_block);
        d = sqrt(d2);
        [weight_sum_per_block(:, pair_idx), weighted_diff2_sum_x_per_block(:, pair_idx), weighted_diff2_sum_y_per_block(:, pair_idx)] = ...
            blockij_variogram_accumulation(d, x(idx_i), x(idx_j), y(idx_i), y(idx_j), spec, is_diagonal_block);
    end
    weight_sum = sum(weight_sum_per_block, 2);
    weighted_diff2_sum_x = sum(weighted_diff2_sum_x_per_block, 2);
    weighted_diff2_sum_y = sum(weighted_diff2_sum_y_per_block, 2);
else
    weight_sum = zeros(spec.M, 1);
    weighted_diff2_sum_x = zeros(spec.M, 1);
    weighted_diff2_sum_y = zeros(spec.M, 1);
    
    for pair_idx = 1:setting.n_unique_block_pairs
        block_i = setting.block_i(pair_idx);
        block_j = setting.block_j(pair_idx);
        idx_i = setting.block_start(block_i):setting.block_end(block_i);
        idx_j = setting.block_start(block_j):setting.block_end(block_j);
        is_diagonal_block = setting.is_diagonal_block(pair_idx);
        d2 = blockij_euclidean_d2(coord(idx_i,:), coord(idx_j,:), coord_norm2(idx_i), coord_norm2(idx_j), is_diagonal_block);
        d = sqrt(d2);
        [weight_sum_i, weighted_diff2_sum_xi, weighted_diff2_sum_yi] = ...
            blockij_variogram_accumulation(d, x(idx_i), x(idx_j), y(idx_i), y(idx_j), spec, is_diagonal_block);
        weight_sum = weight_sum + weight_sum_i;
        weighted_diff2_sum_x = weighted_diff2_sum_x + weighted_diff2_sum_xi;
        weighted_diff2_sum_y = weighted_diff2_sum_y + weighted_diff2_sum_yi;
    end
end
vx = nan(spec.M, 1);
vy = nan(spec.M, 1);
valid_h = weight_sum > 0;
vx(valid_h) = 0.5 * weighted_diff2_sum_x(valid_h) ./ weight_sum(valid_h);
vy(valid_h) = 0.5 * weighted_diff2_sum_y(valid_h) ./ weight_sum(valid_h);
    
end