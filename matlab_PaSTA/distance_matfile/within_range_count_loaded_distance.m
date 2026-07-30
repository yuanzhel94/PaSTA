function count = within_range_count_loaded_distance(distance_matfile, distance_variable, ranges, distance_blocks)
%count = within_range_count_loaded_distance(distance_matfile, distance_variable, range, distance_blocks)
%This function counts the number of data point pairs whose distance within
%the effective ranges, using loaded distance matrix

ranges = ranges(:);
n_ranges = length(ranges);
count = zeros(n_ranges, 1);

if distance_blocks.parallel
    parfor (pair_idx = 1:distance_blocks.n_unique_block_pairs, distance_blocks.n_workers)
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        is_diagonal_block = distance_blocks.is_diagonal_block(pair_idx);
        dij = load_valid_distance_blockij(distance_matfile, ...
                distance_variable, ...
                block_i, ...
                block_j, ...
                distance_blocks);
        count_ij = zeros(n_ranges, 1);
        if ~isempty(dij)
            for range_idx = 1:n_ranges
                count_ij(range_idx) = nnz(dij < ranges(range_idx));
            end
            if ~is_diagonal_block
                count_ij = 2 * count_ij;
            end
        end
        count = count + count_ij;
    end

else
    for pair_idx = 1:distance_blocks.n_unique_block_pairs
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        is_diagonal_block = distance_blocks.is_diagonal_block(pair_idx);
        dij = load_valid_distance_blockij(distance_matfile, ...
                distance_variable, ...
                block_i, ...
                block_j, ...
                distance_blocks);
        count_ij = zeros(n_ranges, 1);
        if ~isempty(dij)
            for range_idx = 1:n_ranges
                count_ij(range_idx) = nnz(dij < ranges(range_idx));
            end
            if ~is_diagonal_block
                count_ij = 2 * count_ij;
            end
        end
        count = count + count_ij;
    end
end

end