function [c_para,b,f,fcov] = fit_stable_fix_exponent(h,v,d_mat,PrecomputeVariance,exponent,nugget)
%This function fit a stable variogram model with predetermined exponent
options = optimoptions('lsqcurvefit', 'Display', 'off');

if isempty(PrecomputeVariance)
    PrecomputeVariance = max(v);
end

f_base = @(b, xdata) b(1) .* (1 - exp(-(xdata./b(2)).^exponent));
f_base_nugget = @(b, xdata) b(1) .* (1 - exp(-(xdata./b(2)).^exponent)) + b(3);
f = @(b, xdata) b(1) .* (1 - exp(-(xdata./b(2)).^b(3))) + b(4);
fcov = @(b,xdata)(xdata>0).*( b(1)  -  b(1) .* (1 - exp(-(xdata./b(2)).^b(3))));

% x0 = [PrecomputeVariance,1];
x0 = [PrecomputeVariance,min(h)];
lb = [0,0];
ub = [2*PrecomputeVariance,inf];

if ~nugget
    b=lsqcurvefit(f_base,x0,h,v,lb,ub,options);
    b = [b,exponent,0];
else
    x0 = [x0,0];
    lb = [lb,0];
    ub = [ub,0.5*PrecomputeVariance];
    b=lsqcurvefit(f_base_nugget,x0,h,v,lb,ub,options);
    b = [b(1:2),exponent,b(3)];
end

c_para=fcov(b,full(d_mat)); 
c_para(logical(eye(size(c_para))))=b(1) + b(4); %set diagonal to variance

end