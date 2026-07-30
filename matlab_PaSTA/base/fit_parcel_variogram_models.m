function bs = fit_parcel_variogram_models(v, parcel_specs, parcel_variance, global_b, setting)
%bs = fit_parcel_variogram_models(v, parcel_specs, parcel_variance, global_b, setting)
%This function fits one stable variogram model for each parcel using the
%global exponent and rescales the sill and nugget to the original variance.

n_parc = numel(parcel_specs);
global_b = global_b(:)';

if n_parc == 1
    bs = global_b;
    return
end

parcel_variance = parcel_variance(:);
bs = zeros(n_parc, 4);

for parc = 1:n_parc
    b = fit_stable_variogram_fixed_exponent(v{parc}, parcel_specs{parc}.h, global_b(3), setting.nugget);
    b(1) = b(1) * parcel_variance(parc);
    b(4) = b(4) * parcel_variance(parc);
    bs(parc, :) = b;
end

end