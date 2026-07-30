function [dmin, dmax] = distance_limits_with_blocks_loaded(distance_matfile, distance_variable, distance_blocks)
%[dmin, dmax] = distance_limits_with_blocks_loaded(distance_matfile, distance_variable, distance_blocks)
%Compute the global min and max distance between points using loaded blocked data

mins_per_block = inf(distance_blocks.n_unique_block_pairs, 1);
maxs_per_block = zeros(distance_blocks.n_unique_block_pairs, 1);

if distance_blocks.parallel
    parfor (pair_idx = 1:distance_blocks.n_unique_block_pairs, distance_blocks.n_workers)
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        
        dij = load_valid_distance_blockij(distance_matfile, ...
            distance_variable, ...
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
    for pair_idx = 1:distance_blocks.n_unique_block_pairs
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        dij = load_valid_distance_blockij(distance_matfile, ...
            distance_variable, ...
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