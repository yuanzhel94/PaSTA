function [dij, idx_i, idx_j] = load_valid_distance_blockij(distance_matfile, distance_variable, block_i, block_j, distance_blocks)
%[dij, idx_i, idx_j] = load_valid_distance_blockij(distance_matfile, distance_variable, block_i, block_j, distance_blocks)
%This function loads the distance block for valid data only

row_start = distance_blocks.block_start(block_i);
row_end = distance_blocks.block_end(block_i);
col_start = distance_blocks.block_start(block_j);
col_end = distance_blocks.block_end(block_j);

dmat_idx_i = (row_start:row_end)';
dmat_idx_j = (col_start:col_end)';

keep_i = distance_blocks.valid(dmat_idx_i);
keep_j = distance_blocks.valid(dmat_idx_j);

if ~any(keep_i) || ~any(keep_j)
    dij = [];
    idx_i = [];
    idx_j = [];
    return
end

idx_i = distance_blocks.in_idx_to_valid_idx(dmat_idx_i(keep_i));
idx_j = distance_blocks.in_idx_to_valid_idx(dmat_idx_j(keep_j));

dij = load_distance_blockij(distance_matfile, ...
    distance_variable, ...
    row_start, row_end, ...
    col_start, col_end);

if ~all(keep_i) || ~all(keep_j)
    dij = dij(keep_i, keep_j);
end

end