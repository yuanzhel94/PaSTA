function count = within_range_count_euclidean(coord, coord_norm2, range, setting)
%count = within_range_count_euclidean(coord, coord_norm2, range, setting)
%This function counts the number of data point pairs whose distance within
%the effective range

range2 = range ^ 2;
if setting.parallel
    count_per_block = zeros(setting.n_unique_block_pairs, 1);
    parfor (pair_idx = 1:setting.n_unique_block_pairs, setting.n_workers)
        block_i = setting.block_i(pair_idx);
        block_j = setting.block_j(pair_idx);
        idx_i = setting.block_start(block_i):setting.block_end(block_i);
        idx_j = setting.block_start(block_j):setting.block_end(block_j);
        is_diagonal_block = setting.is_diagonal_block(pair_idx);
        d2ij = blockij_euclidean_d2(coord(idx_i,:), ...
            coord(idx_j,:), ...
            coord_norm2(idx_i), ...
            coord_norm2(idx_j), ...
            is_diagonal_block);
        count_ij = nnz(d2ij < range2);
        if ~is_diagonal_block
            count_ij = 2 * count_ij;
        end
        count_per_block(pair_idx) = count_ij;
    end
    count = sum(count_per_block);
else
    count = 0;
    for pair_idx = 1:setting.n_unique_block_pairs
        block_i = setting.block_i(pair_idx);
        block_j = setting.block_j(pair_idx);
        idx_i = setting.block_start(block_i):setting.block_end(block_i);
        idx_j = setting.block_start(block_j):setting.block_end(block_j);
        is_diagonal_block = setting.is_diagonal_block(pair_idx);
        d2ij = blockij_euclidean_d2(coord(idx_i,:), ...
            coord(idx_j,:), ...
            coord_norm2(idx_i), ...
            coord_norm2(idx_j), ...
            is_diagonal_block);
        count_ij = nnz(d2ij < range2);
        if ~is_diagonal_block
            count_ij = 2 * count_ij;
        end
        count = count + count_ij;
    end
end
    
end