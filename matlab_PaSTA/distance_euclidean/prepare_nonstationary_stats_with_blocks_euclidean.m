function [rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = prepare_nonstationary_stats_with_blocks_euclidean(coord, coord_norm2, parc_idx_x, bs_x, parc_idx_y, bs_y, dim, setting)
%[rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = prepare_nonstationary_stats_with_blocks_euclidean(coord, coord_norm2, parc_idx_x, bs_x, parc_idx_y, bs_y, dim, setting)
%Prepare stats from each block's distance matrix and nonstationary covariance model, for
%future nef computation.

trace_x = sum(bs_x(parc_idx_x, 1) + bs_x(parc_idx_x, 4));
trace_y = sum(bs_y(parc_idx_y, 1) + bs_y(parc_idx_y, 4));

rowsum_x = zeros(setting.N, 1);
rowsum_y = zeros(setting.N, 1);
cov_prodsum = 0;

if setting.parallel
    parfor (block_i = 1:setting.n_blocks, setting.n_workers)
        rowsum_xi = zeros(setting.N, 1);
        rowsum_yi = zeros(setting.N, 1);
        cov_prodsumi = 0;
        idx_i = setting.block_start(block_i):setting.block_end(block_i);

        bxi = bs_x(parc_idx_x(idx_i), :);
        byi = bs_y(parc_idx_y(idx_i), :);

        for block_j = block_i:setting.n_blocks
            idx_j = setting.block_start(block_j):setting.block_end(block_j);
            is_diagonal_block = block_i == block_j;
            d2ij = blockij_euclidean_d2(coord(idx_i, :), ...
                coord(idx_j, :), ...
                coord_norm2(idx_i), ...
                coord_norm2(idx_j), ...
                is_diagonal_block);
            dij = sqrt(d2ij);

            bxj = bs_x(parc_idx_x(idx_j), :);
            byj = bs_y(parc_idx_y(idx_j), :);

            covij_x = blockij_nonstationary_covariance(dij, bxi, bxj, dim, is_diagonal_block);
            covij_y = blockij_nonstationary_covariance(dij, byi, byj, dim, is_diagonal_block);
            
            rowsum_xi(idx_i) = rowsum_xi(idx_i) + sum(covij_x, 2);
            rowsum_yi(idx_i) = rowsum_yi(idx_i) + sum(covij_y, 2);
            if is_diagonal_block
                cov_prodsumi = cov_prodsumi + dot(covij_x(:), covij_y(:));
            else
                rowsum_xi(idx_j) = rowsum_xi(idx_j) + sum(covij_x, 1)';
                rowsum_yi(idx_j) = rowsum_yi(idx_j) + sum(covij_y, 1)';
                cov_prodsumi = cov_prodsumi + 2 * dot(covij_x(:), covij_y(:));
            end
        end
        rowsum_x = rowsum_x + rowsum_xi;
        rowsum_y = rowsum_y + rowsum_yi;
        cov_prodsum = cov_prodsum + cov_prodsumi;
    end
else
    for block_i = 1:setting.n_blocks
        idx_i = setting.block_start(block_i):setting.block_end(block_i);

        bxi = bs_x(parc_idx_x(idx_i), :);
        byi = bs_y(parc_idx_y(idx_i), :);
        
        for block_j = block_i:setting.n_blocks
            idx_j = setting.block_start(block_j):setting.block_end(block_j);
            is_diagonal_block = block_i == block_j;
            d2ij = blockij_euclidean_d2(coord(idx_i, :), ...
                coord(idx_j, :), ...
                coord_norm2(idx_i), ...
                coord_norm2(idx_j), ...
                is_diagonal_block);
            dij = sqrt(d2ij);

            bxj = bs_x(parc_idx_x(idx_j), :);
            byj = bs_y(parc_idx_y(idx_j), :);

            covij_x = blockij_nonstationary_covariance(dij, bxi, bxj, dim, is_diagonal_block);
            covij_y = blockij_nonstationary_covariance(dij, byi, byj, dim, is_diagonal_block);
            
            rowsum_x(idx_i) = rowsum_x(idx_i) + sum(covij_x, 2);
            rowsum_y(idx_i) = rowsum_y(idx_i) + sum(covij_y, 2);
            if is_diagonal_block
                cov_prodsum = cov_prodsum + dot(covij_x(:), covij_y(:));
            else
                rowsum_x(idx_j) = rowsum_x(idx_j) + sum(covij_x, 1)';
                rowsum_y(idx_j) = rowsum_y(idx_j) + sum(covij_y, 1)';
                cov_prodsum = cov_prodsum + 2 * dot(covij_x(:), covij_y(:));
            end
        end
    end
end


end