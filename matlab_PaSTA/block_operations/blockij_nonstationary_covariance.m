function covij = blockij_nonstationary_covariance(dij, bi, bj, dim, is_diagonal_block)
%covij = blockij_nonstationary_covariance(d, bi, bj, dim, is_diagonal_block)
%This funciton compute the nonstationary covariance entries for block ij

sill_i = bi(:,1);
sill_j = bj(:,1);
phii = bi(:,2);
phij = bj(:,2);
exponent = bi(1, 3); % exponent is global

sig = (phii .^ 2 + (phij .^ 2)') ./ 2; % (ni, nj)
qij = dij .^ 2 ./ sig; % (ni, nj)
covij = ((phii * phij') ./ sig) .^ (dim / 2) .* (sqrt(sill_i) * sqrt(sill_j')) .* exp(- (sqrt(qij) .^ exponent) );

if is_diagonal_block
    n = size(covij, 1);
    covij(1:(n + 1):end) = bi(:, 1) + bi(:, 4);
end

end