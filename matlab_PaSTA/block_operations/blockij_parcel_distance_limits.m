function [dmin, dmax] = blockij_parcel_distance_limits(dij, parc_i, parc_j, n_parc)
%[dmin, dmax] = blockij_parcel_distance_limits(dij, parc_i, parc_j, n_parc)
%This function computes within parcel min and max distance for one loaded
%distance block

dmin = inf(n_parc, 1);
dmax = zeros(n_parc, 1);
for parc = 1:n_parc
    selected_i = parc_i == parc;
    selected_j = parc_j == parc;
    if ~any(selected_i) || ~any(selected_j)
        continue
    end
    dij_parc = dij(selected_i, selected_j);
    dmax(parc) = max(dij_parc, [], 'all');
    dij_parc(dij_parc <= 0) = Inf;
    dmin(parc) = min(dij_parc, [], 'all');
end

end