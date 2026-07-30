function [rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = prepare_nonstationary_stats_with_blocks_in_memory_distance(distance_memory, parc_idx_x, bs_x, parc_idx_y, bs_y, dim, distance_blocks)
%[rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = prepare_nonstationary_stats_with_blocks_in_memory_distance(distance_memory, parc_idx_x, bs_x, parc_idx_y, bs_y, dim, distance_blocks)
%This function prepares covariance statistics using an in-memory full
%distance matrix or strict upper-triangular distance vector and
%nonstationary parcel-level covariance models in blockwise manner.

parc_idx_x = parc_idx_x(:);
parc_idx_y = parc_idx_y(:);
N = nnz(distance_blocks.valid);

trace_x = sum(bs_x(parc_idx_x, 1) + bs_x(parc_idx_x, 4));
trace_y = sum(bs_y(parc_idx_y, 1) + bs_y(parc_idx_y, 4));

rowsum_x = zeros(N, 1);
rowsum_y = zeros(N, 1);
cov_prodsum = 0;

if distance_blocks.parallel
    parfor (block_i = 1:distance_blocks.n_blocks, distance_blocks.n_workers)
        rowsum_xi = zeros(N, 1);
        rowsum_yi = zeros(N, 1);
        cov_prodsumi = 0;

        for block_j = block_i:distance_blocks.n_blocks
            is_diagonal_block = block_i == block_j;
            [dij, idx_i, idx_j] = load_valid_distance_blockij_in_memory(distance_memory, block_i, block_j, distance_blocks);
            if isempty(dij)
                continue
            end

            bxi = bs_x(parc_idx_x(idx_i), :);
            bxj = bs_x(parc_idx_x(idx_j), :);
            byi = bs_y(parc_idx_y(idx_i), :);
            byj = bs_y(parc_idx_y(idx_j), :);

            covij_x = blockij_nonstationary_covariance(dij, bxi, bxj, dim, is_diagonal_block);
            covij_y = blockij_nonstationary_covariance(dij, byi, byj, dim, is_diagonal_block);

            rowsum_xi(idx_i) = rowsum_xi(idx_i) + sum(covij_x, 2);
            rowsum_yi(idx_i) = rowsum_yi(idx_i) + sum(covij_y, 2);

            if is_diagonal_block
                cov_prodsumi = cov_prodsumi + dot(covij_x(:), covij_y(:));
            else
                rowsum_xi(idx_j) = rowsum_xi(idx_j) + sum(covij_x, 1).';
                rowsum_yi(idx_j) = rowsum_yi(idx_j) + sum(covij_y, 1).';
                cov_prodsumi = cov_prodsumi + 2 * dot(covij_x(:), covij_y(:));
            end
        end

        rowsum_x = rowsum_x + rowsum_xi;
        rowsum_y = rowsum_y + rowsum_yi;
        cov_prodsum = cov_prodsum + cov_prodsumi;
    end

else
    for block_i = 1:distance_blocks.n_blocks
        for block_j = block_i:distance_blocks.n_blocks
            is_diagonal_block = block_i == block_j;

            [dij, idx_i, idx_j] = load_valid_distance_blockij_in_memory(distance_memory, block_i, block_j, distance_blocks);

            if isempty(dij)
                continue
            end

            bxi = bs_x(parc_idx_x(idx_i), :);
            bxj = bs_x(parc_idx_x(idx_j), :);
            byi = bs_y(parc_idx_y(idx_i), :);
            byj = bs_y(parc_idx_y(idx_j), :);

            covij_x = blockij_nonstationary_covariance(dij, bxi, bxj, dim, is_diagonal_block);
            covij_y = blockij_nonstationary_covariance(dij, byi, byj, dim, is_diagonal_block);

            rowsum_x(idx_i) = rowsum_x(idx_i) + sum(covij_x, 2);
            rowsum_y(idx_i) = rowsum_y(idx_i) + sum(covij_y, 2);

            if is_diagonal_block
                cov_prodsum = cov_prodsum + dot(covij_x(:), covij_y(:));
            else
                rowsum_x(idx_j) = rowsum_x(idx_j) + sum(covij_x, 1).';
                rowsum_y(idx_j) = rowsum_y(idx_j) + sum(covij_y, 1).';
                cov_prodsum = cov_prodsum + 2 * dot(covij_x(:), covij_y(:));
            end
        end
    end
end

end