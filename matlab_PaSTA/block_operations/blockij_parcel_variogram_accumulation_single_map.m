function [weight_sum, weighted_diff2_sum] = blockij_parcel_variogram_accumulation_single_map(dij, xi, xj, parc_i, parc_j, parcel_specs, is_diagonal_block)
%[weight_sum, weighted_diff2_sum] = blockij_parcel_variogram_accumulation_single_map(dij, xi, xj, parc_i, parc_j, parcel_specs, is_diagonal_block)
%This function prepares data for variogram computation at parcel level for
%one map using loaded distance block ij.

n_parc = numel(parcel_specs);

weight_sum = cell(n_parc, 1);
weighted_diff2_sum = cell(n_parc, 1);

for parc = 1:n_parc
    spec = parcel_specs{parc};
    weight_sum{parc} = zeros(spec.M, 1);
    weighted_diff2_sum{parc} = zeros(spec.M, 1);
    selected_i = parc_i == parc;
    selected_j = parc_j == parc;
    if ~any(selected_i) || ~any(selected_j)
        continue
    end
    [weight_sum{parc}, weighted_diff2_sum{parc}] = blockij_variogram_accumulation_single_map(dij(selected_i, selected_j), ...
            xi(selected_i), ...
            xj(selected_j), ...
            spec, ...
            is_diagonal_block);
end

end