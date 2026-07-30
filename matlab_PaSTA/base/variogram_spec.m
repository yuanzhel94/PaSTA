function spec = variogram_spec(dmin, dmax, M, qd, kernel_scale, nugget)
%spec = variogram_spec(dmin, dmax, M, qd)
%This function derive a struct that contains the specifications for
%variogram estimation.

variogram_dmax = qd * dmax;
if variogram_dmax <= dmin
    error('qd*dmax must be greater than dmin.');
end

spec = struct();
spec.M = double(M);
spec.qd = double(qd);
spec.kernel_scale = double(kernel_scale);
spec.nugget = logical(nugget);

spec.dmin = double(dmin);
spec.dmax = double(dmax);
spec.variogram_dmin = double(dmin);
spec.variogram_dmax = double(variogram_dmax);

spec.h = linspace(spec.variogram_dmin, spec.variogram_dmax, spec.M)';
spec.lag_sep = (spec.variogram_dmax - spec.variogram_dmin) / (M - 1); 
spec.delta = spec.lag_sep / 2;
spec.sigma = 6 * spec.delta;

spec.L = kernel_scale * spec.sigma / 2.68; % kernel (half) size to truncate weight computation, i.e., ignore pairs whose distance > 4*sigma/2.68 from bin center.
spec.L2 = spec.L ^ 2;
spec.weight_constant = - 2.68 ^ 2 / (2 * spec.sigma ^ 2); % constant term in kernel weight before exponential

spec.kernel_nh = ceil(spec.L / spec.lag_sep);

spec.bin_idx_offset = -spec.kernel_nh:spec.kernel_nh; % offset idx of neighbouring bins fail within L

end