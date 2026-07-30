function d2 = blockij_euclidean_d2(coord_i, coord_j, coord_i_norm2, coord_j_norm2, is_diagonal_block)
% d2 = blockij_euclidean_d2(coord_i, coord_j, coord_i_norm2, coord_j_norm2, is_diagonal_block)
% This function computes the distance between points in two blocks
% from Inputs:
%   coord_i: (Ni, dim) euclidean space coordinates for block i
%   coord_j: (Nj, dim) euclidean space coordinates for block j
%   coord_i_norm2: (Ni, 1) squared L2 norm of coord_i, precompute for computational efficiency when loop over
%   coord_j_norm2: (Nj, 1) squared L2 norm of coord_j, precompute for computational efficiency when loop over
%   is_diagonal_block: bool indicator of if (i,j) is diagonal block, i.e., i==j

d2 = coord_i * coord_j'; % (Ni,Nj) product
d2 = -2 .* d2;
d2 = d2 + coord_i_norm2;
d2 = d2 + coord_j_norm2'; % final squared distance
d2 = max(d2, 0); % avoid floating point inaccuracy for small separation - unrealistic edge case for realworld neuroimaging application
if is_diagonal_block
    n = size(d2, 1);
    d2(1:(n + 1):end) = 0; % avoid floating point inaccuracy for self distance
end

end