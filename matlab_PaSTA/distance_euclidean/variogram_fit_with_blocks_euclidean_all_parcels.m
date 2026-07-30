function bs = variogram_fit_with_blocks_euclidean_all_parcels(x, coord, coord_norm2, parc_idx, n_parc, global_b, setting)
%bs = variogram_fit_with_blocks_euclidean_all_parcels(x, coord, coord_norm2, parc_idx, global_b, setting)
%This function estimate and fit the variogram for derived parcels, in a
%blockwise manner for memory efficiency

if n_parc == 1
    bs = global_b;
    return
end

bs = zeros(n_parc, 4);

if setting.parallel
    parfor (parc_i = 1:n_parc, setting.n_workers)
        selected = parc_idx == parc_i;
        bs(parc_i, :) = variogram_fit_with_blocks_euclidean_single_parcel(x(selected), coord(selected,:), coord_norm2(selected), global_b, setting);
    end
else
    for parc_i = 1:n_parc
        selected = parc_idx == parc_i;
        bs(parc_i, :) = variogram_fit_with_blocks_euclidean_single_parcel(x(selected), coord(selected,:), coord_norm2(selected), global_b, setting);
    end
end

end