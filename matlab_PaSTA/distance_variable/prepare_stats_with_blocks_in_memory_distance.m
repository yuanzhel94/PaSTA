function [rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = prepare_stats_with_blocks_in_memory_distance(distance_memory, bx, by, distance_blocks)
%[rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum] = prepare_stats_with_blocks_in_memory_distance(distance_memory, bx, by, distance_blocks)
%This function prepares covariance statistics using an in-memory full
%distance matrix or strict upper-triangular distance vector in blockwise
%manner.

N = nnz(distance_blocks.valid);

trace_x = N * (bx(1) + bx(4));
trace_y = N * (by(1) + by(4));
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
            covij_x = blockij_stable_covariance(dij, bx, is_diagonal_block);
            covij_y = blockij_stable_covariance(dij, by, is_diagonal_block);
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
            covij_x = blockij_stable_covariance(dij, bx, is_diagonal_block);
            covij_y = blockij_stable_covariance(dij, by, is_diagonal_block);
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