function [x_standardized, parcel_variance, parcel_specs] = prepare_parcel_variogram_inputs(x, parc_idx, n_parc, dmin, dmax, setting)
%[x_standardized, parcel_variance, parcel_specs] = prepare_parcel_variogram_inputs(x, parc_idx, n_parc, dmin, dmax, setting)
%This function standardizes one map within each parcel and creates the
%parcel-specific variogram specifications.

x = x(:);
parc_idx = parc_idx(:);
dmin = dmin(:);
dmax = dmax(:);

x_standardized = zeros(size(x));
parcel_variance = zeros(n_parc, 1);
parcel_specs = cell(n_parc, 1);

for parc = 1:n_parc
    selected = parc_idx == parc;
    n_parcel = nnz(selected);
    x_parcel = x(selected);
    variance_parcel = var(x_parcel);
    if ~isfinite(variance_parcel) || variance_parcel <= 0
        error('Parcel %d has zero or invalid variance.', parc);
    end
    x_standardized(selected) = (x_parcel - mean(x_parcel)) ./ sqrt(variance_parcel);
    parcel_variance(parc) = variance_parcel;
    M_parcel = setting.M_scale * ceil(sqrt(n_parcel));
    parcel_specs{parc} = variogram_spec(dmin(parc), ...
        dmax(parc), ...
        M_parcel, ...
        setting.qd, ...
        setting.kernel_scale, ...
        setting.nugget);
end

end