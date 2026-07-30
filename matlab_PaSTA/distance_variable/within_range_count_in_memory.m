function count = within_range_count_in_memory(ranges, distance_memory, distance_blocks)
%count = within_range_count_in_memory(ranges, distance_memory, distance_blocks)
%This function counts the number of distance matrix entries less than or
%equal to each range using blockwise traversal of an in-memory full distance
%matrix or strict upper-triangular distance vector.

ranges = ranges(:);
n_ranges = numel(ranges);
count = zeros(n_ranges, 1);

if distance_blocks.parallel
    parfor (pair_idx = 1:distance_blocks.n_unique_block_pairs, distance_blocks.n_workers)
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        dij = load_valid_distance_blockij_in_memory(distance_memory, block_i, block_j, distance_blocks);
        count_ij = zeros(n_ranges, 1);
        if ~isempty(dij)
            for range_idx = 1:n_ranges
                if block_i == block_j
                    count_ij(range_idx) = nnz(dij <= ranges(range_idx));
                else
                    count_ij(range_idx) = 2 * nnz(dij <= ranges(range_idx));
                end
            end
        end
        count = count + count_ij;
    end

else
    for pair_idx = 1:distance_blocks.n_unique_block_pairs
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        dij = load_valid_distance_blockij_in_memory(distance_memory, block_i, block_j, distance_blocks);
        count_ij = zeros(n_ranges, 1);
        if ~isempty(dij)
            for range_idx = 1:n_ranges
                if block_i == block_j
                    count_ij(range_idx) = nnz(dij <= ranges(range_idx));
                else
                    count_ij(range_idx) = 2 * nnz(dij <= ranges(range_idx));
                end
            end
        end
        count = count + count_ij;
    end
end

end