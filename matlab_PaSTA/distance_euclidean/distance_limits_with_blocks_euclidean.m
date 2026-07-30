function [dmin, dmax] = distance_limits_with_blocks_euclidean(coord, coord_norm2, setting)
%[dmin, dmax] = distance_limits_with_blocks(coord, coord_norm2, setting)
%Compute the global min and max distance between points using blocked data

mins_d2_per_block = inf(setting.n_unique_block_pairs, 1);
maxs_d2_per_block = zeros(setting.n_unique_block_pairs, 1);

if setting.parallel
    parfor (pair_idx = 1:setting.n_unique_block_pairs, setting.n_workers)
        block_i = setting.block_i(pair_idx);
        block_j = setting.block_j(pair_idx);
        idx_i = setting.block_start(block_i):setting.block_end(block_i);
        idx_j = setting.block_start(block_j):setting.block_end(block_j);
        [mins_d2_per_block(pair_idx), maxs_d2_per_block(pair_idx)] = blockij_euclidean_d2_limits(...
            coord, ...
            coord_norm2, ...
            idx_i, ...
            idx_j, ...
            setting.is_diagonal_block(pair_idx));
    end
else
    for pair_idx = 1:setting.n_unique_block_pairs
        block_i = setting.block_i(pair_idx);
        block_j = setting.block_j(pair_idx);
        idx_i = setting.block_start(block_i):setting.block_end(block_i);
        idx_j = setting.block_start(block_j):setting.block_end(block_j);
        [mins_d2_per_block(pair_idx), maxs_d2_per_block(pair_idx)] = blockij_euclidean_d2_limits(...
            coord, ...
            coord_norm2, ...
            idx_i, ...
            idx_j, ...
            setting.is_diagonal_block(pair_idx));
    end
end
d2min = min(mins_d2_per_block);
d2max = max(maxs_d2_per_block);
dmin = sqrt(d2min);
dmax = sqrt(d2max);
end