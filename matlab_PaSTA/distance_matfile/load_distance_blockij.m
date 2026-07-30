function dij = load_distance_blockij(distance_matfile, distance_variable, row_start, row_end, col_start, col_end)
%dij = load_distance_blockij(distance_matfile, distance_variable, row_start, row_end, col_start, col_end)
%This function loads partial distance matrix from precomputed matfile,
%matfile needs to be v7.3 to support partial loading of huge matrix
%

row_idx = row_start:row_end;
col_idx = col_start:col_end;
file_idx = substruct('.', distance_variable, '()', {row_idx, col_idx});
dij = double(subsref(distance_matfile, file_idx));

end