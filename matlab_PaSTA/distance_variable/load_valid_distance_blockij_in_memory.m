function [dij, idx_i, idx_j] = load_valid_distance_blockij_in_memory(distance_memory, block_i, block_j, distance_blocks)
%[dij, idx_i, idx_j] = load_valid_distance_blockij_in_memory(distance_memory, block_i, block_j, distance_blocks)
%This function extracts one in-memory distance block, removes observations
%that are invalid in x, y, or coordinates, and returns indices in the
%filtered valid-data space.
%
%dij:
%   Double-precision distance block containing valid observations only.
%
%idx_i and idx_j:
%   Indices of the rows and columns of dij in the filtered valid-data space.

row_start = distance_blocks.block_start(block_i);
row_end = distance_blocks.block_end(block_i);
col_start = distance_blocks.block_start(block_j);
col_end = distance_blocks.block_end(block_j);

distance_idx_i = (row_start:row_end)';
distance_idx_j = (col_start:col_end)';

keep_i = distance_blocks.valid(distance_idx_i);
keep_j = distance_blocks.valid(distance_idx_j);

if ~any(keep_i) || ~any(keep_j)
    dij = [];
    idx_i = [];
    idx_j = [];
    return
end

idx_i = distance_blocks.in_idx_to_valid_idx(distance_idx_i(keep_i));
idx_j = distance_blocks.in_idx_to_valid_idx(distance_idx_j(keep_j));

dij = load_distance_blockij_in_memory( ...
    distance_memory, ...
    block_i, ...
    block_j, ...
    distance_blocks);

if ~all(keep_i) || ~all(keep_j)
    dij = dij(keep_i, keep_j);
end

end