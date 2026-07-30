function [vx, vy] = global_variogram_with_blocks_in_memory(x, y, distance_memory, distance_blocks, spec)
%[vx, vy] = global_variogram_with_blocks_in_memory(x, y, distance_memory, distance_blocks, spec)
%This function computes global empirical variograms using blockwise
%traversal of an in-memory full distance matrix or strict upper-triangular
%distance vector.

weight_sum = zeros(spec.M, 1);
weighted_diff2_sum_x = zeros(spec.M, 1);
weighted_diff2_sum_y = zeros(spec.M, 1);

if distance_blocks.parallel
    parfor (pair_idx = 1:distance_blocks.n_unique_block_pairs, distance_blocks.n_workers)
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        [dij, idx_i, idx_j] = load_valid_distance_blockij_in_memory(distance_memory, ...
                block_i, ...
                block_j, ...
                distance_blocks);
        weight_sum_ij = zeros(spec.M, 1);
        weighted_diff2_sum_xij = zeros(spec.M, 1);
        weighted_diff2_sum_yij = zeros(spec.M, 1);
        if ~isempty(dij)
            is_diagonal_block = block_i == block_j;
            [weight_sum_ij, weighted_diff2_sum_xij, weighted_diff2_sum_yij] = ...
                blockij_variogram_accumulation(dij, ...
                    x(idx_i), ...
                    x(idx_j), ...
                    y(idx_i), ...
                    y(idx_j), ...
                    spec, ...
                    is_diagonal_block);
        end

        weight_sum = weight_sum + weight_sum_ij;
        weighted_diff2_sum_x = weighted_diff2_sum_x + weighted_diff2_sum_xij;
        weighted_diff2_sum_y = weighted_diff2_sum_y + weighted_diff2_sum_yij;
    end

else
    for pair_idx = 1:distance_blocks.n_unique_block_pairs
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        [dij, idx_i, idx_j] = load_valid_distance_blockij_in_memory(distance_memory, ...
                block_i, ...
                block_j, ...
                distance_blocks);
        weight_sum_ij = zeros(spec.M, 1);
        weighted_diff2_sum_xij = zeros(spec.M, 1);
        weighted_diff2_sum_yij = zeros(spec.M, 1);
        if ~isempty(dij)
            is_diagonal_block = block_i == block_j;
            [weight_sum_ij, weighted_diff2_sum_xij, weighted_diff2_sum_yij] = ...
                blockij_variogram_accumulation(dij, ...
                    x(idx_i), ...
                    x(idx_j), ...
                    y(idx_i), ...
                    y(idx_j), ...
                    spec, ...
                    is_diagonal_block);
        end

        weight_sum = weight_sum + weight_sum_ij;
        weighted_diff2_sum_x = weighted_diff2_sum_x + weighted_diff2_sum_xij;
        weighted_diff2_sum_y = weighted_diff2_sum_y + weighted_diff2_sum_yij;
    end
end

vx = nan(spec.M, 1);
vy = nan(spec.M, 1);

valid_bin = weight_sum > 0;

vx(valid_bin) = 0.5 .* weighted_diff2_sum_x(valid_bin) ./ weight_sum(valid_bin);

vy(valid_bin) = 0.5 .* weighted_diff2_sum_y(valid_bin) ./ weight_sum(valid_bin);

end