function [v,h] = estimate_variogram(D, data, M, qd)
% Estimate the empirical variogram from distance matrix between vertices,
% and data value at each vertex. Estiamtion performed in M bins, ranging
% from min_distance to qd * max_distance, where max_distance is the max distance in
% the distance matrix.
%
%D:
% Distance matrix between all vertices
%data:
% Data value at each vertex.
%M:
% number of bins to estimate variogram.
%qd:
% maximum distance to evaluate variogram
%
% Returns:
% v: semivariance of the varigoram (M,1)
% h: lag distances of the variogram (1, M)

Dmax = qd * max(D(:));
Dmin = min(D(D > 0));

mask = triu(D > 0 & D <= Dmax, 1);

[row, col] = find(mask);
dval = D(mask);

diff2 = (data(row) - data(col)).^2;

% Sort by distance
[dval, order] = sort(dval);
diff2 = diff2(order);

h = linspace(Dmin, Dmax, M);
delta = (Dmax - Dmin) / (M - 1) * 0.5;
sigma = 6 * delta;

L = 4 * sigma / 2.68;

v = zeros(M, 1);

for i = 1:M
    lo = h(i) - L;
    hi = h(i) + L;

    left = lower_bound(dval, lo);
    right = upper_bound(dval, hi);

    if right >= left
        idx = left:right;

        local_d = dval(idx);
        local_diff2 = diff2(idx);

        w = exp(-((2.68 * abs(h(i) - local_d)).^2) / (2 * sigma^2));

        v(i) = 0.5 * sum(w .* local_diff2) / sum(w);
    else
        v(i) = NaN;
    end
end

end
