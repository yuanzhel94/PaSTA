function [dmin, dmax] = distance_limits_with_blocks_in_memory(distance_memory, distance_blocks)
%[dmin, dmax] = distance_limits_with_blocks_in_memory(distance_memory, distance_blocks)
%This function computes the global minimum positive distance and maximum
%distance using blockwise traversal of an in-memory distance matrix or
%strict upper-triangular distance vector.

n_pairs = distance_blocks.n_unique_block_pairs;

mins_per_block = inf(n_pairs, 1);
maxs_per_block = zeros(n_pairs, 1);

if distance_blocks.parallel
    parfor (pair_idx = 1:n_pairs, distance_blocks.n_workers)
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        dij = load_valid_distance_blockij_in_memory(distance_memory, ...
            block_i, ...
            block_j, ...
            distance_blocks);
        if isempty(dij)
            continue
        end
        maxs_per_block(pair_idx) = max(dij, [], 'all');
        dij(dij <= 0) = Inf;
        mins_per_block(pair_idx) = min(dij, [], 'all');
    end

else
    for pair_idx = 1:n_pairs
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        dij = load_valid_distance_blockij_in_memory(distance_memory, ...
            block_i, ...
            block_j, ...
            distance_blocks);
        if isempty(dij)
            continue
        end
        maxs_per_block(pair_idx) = max(dij, [], 'all');
        dij(dij <= 0) = Inf;
        mins_per_block(pair_idx) = min(dij, [], 'all');
    end
end

dmin = min(mins_per_block);
dmax = max(maxs_per_block);

end