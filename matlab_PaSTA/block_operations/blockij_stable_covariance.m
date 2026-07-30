function covij = blockij_stable_covariance(d, b, is_diagonal_block)
%covij = blockij_stable_covariance(d, b, is_diagonal_block)
%This funciton compute the covariance entries for block ij

covij = b(1) .* exp(-(d ./ b(2)) .^ b(3));
if is_diagonal_block
    n = size(covij, 1);
    covij(1:(n+1):end) = b(1) + b(4); % set diagonal to sill + nugget
end

end