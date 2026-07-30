function [b, f] = fit_stable_variogram_fixed_exponent(v, h, exponent, nugget)
%[b, f] = fit_stable_variogram_fixed_exponent(v, h, exponent, nugget)
%%This function fit a stable variogram model with predetermined exponent

valid = isfinite(h) & isfinite(v) & h > 0 & v >= 0;
h_fit = h(valid);
v_fit = v(valid);
h_min = min(h_fit);
v_max = max(v_fit);

options = optimoptions('lsqcurvefit', 'Display', 'off');

f = @(b, xdata) b(1) .* (1 - exp(-(xdata ./ b(2)) .^ b(3))) + b(4);
x0 = [v_max,h_min];
lb = [0,0];
ub = [2*v_max,inf];
if ~nugget
    f_fit = @(b, xdata) b(1) .* (1 - exp(-(xdata ./ b(2)) .^ exponent));
    b=lsqcurvefit(f_fit,x0,h_fit,v_fit,lb,ub,options);
    b = [b,exponent,0];
else
    f_fit = @(b, xdata) b(1) .* (1 - exp(-(xdata ./ b(2)) .^ exponent)) + b(3);
    x0 = [x0,0];
    lb = [lb,0];
    ub = [ub,0.5*v_max];
    b=lsqcurvefit(f_fit,x0,h_fit,v_fit,lb,ub,options);
    b = [b(1:2),exponent,b(3)];
end

end