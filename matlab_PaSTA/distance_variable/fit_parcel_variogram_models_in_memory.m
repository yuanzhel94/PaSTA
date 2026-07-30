function [bs_x, bs_y] = fit_parcel_variogram_models_in_memory(x, y, distance_memory, parc_idx_x, n_parc_x, parc_idx_y, n_parc_y, bx, by, setting, distance_blocks)
%[bs_x, bs_y] = fit_parcel_variogram_models_in_memory(x, y, distance_memory, parc_idx_x, n_parc_x, parc_idx_y, n_parc_y, bx, by, setting, distance_blocks)
%This function fits parcel-level stable variogram models for two maps using
%an in-memory full distance matrix or strict upper-triangular distance
%vector in blockwise manner

x = x(:);
y = y(:);
parc_idx_x = parc_idx_x(:);
parc_idx_y = parc_idx_y(:);

if n_parc_x == 1 && n_parc_y == 1
    bs_x = bx(:)';
    bs_y = by(:)';
    return
end

[dmin_x, dmax_x, dmin_y, dmax_y] = parcel_distance_limits_with_blocks_in_memory(parc_idx_x, n_parc_x, parc_idx_y, n_parc_y, distance_memory, distance_blocks);

[x_standardized, parcel_variance_x, parcel_specs_x] = prepare_parcel_variogram_inputs(x, parc_idx_x, n_parc_x, dmin_x, dmax_x, setting);

[y_standardized, parcel_variance_y, parcel_specs_y] = prepare_parcel_variogram_inputs(y, parc_idx_y, n_parc_y, dmin_y, dmax_y, setting);

[vx, vy] = parcel_variogram_with_blocks_in_memory(x_standardized, y_standardized, parc_idx_x, parc_idx_y, parcel_specs_x, parcel_specs_y, distance_memory, distance_blocks);

bs_x = fit_parcel_variogram_models(vx, parcel_specs_x, parcel_variance_x, bx, setting);
bs_y = fit_parcel_variogram_models(vy, parcel_specs_y, parcel_variance_y, by, setting);

end