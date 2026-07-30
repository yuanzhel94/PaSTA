function [b, f] = fit_stable_variogram(v, h, nugget)
%[b, f] = fit_stable_variogram(v, h, nugget)
%This function fits the stable variogram model for maps

valid = isfinite(h) & isfinite(v) & h > 0 & v >= 0;
h_fit = h(valid);
v_fit = v(valid);
h_min = min(h_fit);
v_max = max(v_fit);

options = optimoptions('lsqcurvefit', 'Display', 'off');

f_no_nugget = @(b, xdata) b(1) .* (1 - exp(-(xdata ./ b(2)) .^ b(3)));
f = @(b, xdata) b(1) .* (1 - exp(-(xdata ./ b(2)) .^ b(3))) + b(4);

x0 = [v_max,h_min,1];
lb = [0,0,0];
ub = [2*v_max,inf,2];

if ~nugget
    b=lsqcurvefit(f_no_nugget,x0,h_fit,v_fit,lb,ub,options);
    b = [b,0];
else
    x0 = [x0,0];
    lb = [lb,0];
    ub = [ub,0.5*v_max];
    b=lsqcurvefit(f,x0,h_fit,v_fit,lb,ub,options);
end

end