function dij = load_distance_blockij_in_memory(distance_memory, block_i, block_j, distance_blocks)
%dij = load_distance_blockij_in_memory(distance_memory, block_i, block_j, distance_blocks)
%This function extracts one distance block from an in-memory full matrix or
%strict upper-triangular vector. The returned block is double precision.
%
%Only upper-triangular block pairs, where block_j >= block_i, are expected.

if block_j < block_i
    error('block_j must be greater than or equal to block_i.');
end

row_start = distance_blocks.block_start(block_i);
row_end = distance_blocks.block_end(block_i);
col_start = distance_blocks.block_start(block_j);
col_end = distance_blocks.block_end(block_j);

n_rows = row_end - row_start + 1;
n_cols = col_end - col_start + 1;

is_diagonal_block = block_i == block_j;

switch distance_memory.storage
    case 'matrix'
        dij = double(distance_memory.data( ...
            row_start:row_end, ...
            col_start:col_end));
    case 'upper'
        dij = zeros(n_rows, n_cols);
        if is_diagonal_block
            %Reconstruct a complete symmetric diagonal block from the
            %packed strict upper-triangular vector.
            for local_j = 2:n_cols
                global_j = col_start + local_j - 1;
                global_i = (row_start:(global_j - 1))';
                packed_idx = (global_j - 1) * (global_j - 2) / 2 + global_i;
                values = double(distance_memory.data(packed_idx));
                values = values(:);
                dij(1:(local_j - 1), local_j) = values;
                dij(local_j, 1:(local_j - 1)) = values.';
            end

        else
            %For an upper off-diagonal block, every global row index is
            %smaller than every global column index.
            global_i = (row_start:row_end)';
            for local_j = 1:n_cols
                global_j = col_start + local_j - 1;
                packed_idx = (global_j - 1) * (global_j - 2) / 2 + global_i;
                values = double(distance_memory.data(packed_idx));
                dij(:, local_j) = values(:);
            end
        end

    otherwise
        error('Unsupported in-memory distance storage type: %s', ...
            distance_memory.storage);
end

end