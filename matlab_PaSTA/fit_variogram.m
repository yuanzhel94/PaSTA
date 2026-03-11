function [c_para,b,f,fcov]=fit_variogram(h,v,d_mat,fitype,PrecomputeVariance,nugget)
%Fit stable or Matern variogram function to emprical variogram using 
%curve fitting. Return function handles to the fitted variogram (f) and the
%covariance function (fcov). 
%b is a vector of fitted parameters
%c_para is the fitted covariance matrix 
%%%
% h,v:        
%  Empirical variogram. v(h) is the variogram value at a distance of h. 
% d_mat:
%  Distance matrix. Distance between points on the grid. 
% fitype:
%  'stable' or 'matern' 
% PrecomputeVariance
%   If empty, the sill of the fitted variogram is fitted; otherwise the
%   sill is set to the value of PrecomputeVariance. The sill is value of
%   the plateau in the variogram. PrecomputeVariance could be set to the
%   variance across all elements of the map, although this may be a poor
%   approximation in cases of strong autocorrelation, in which case the
%   sill should be fitted. 
%
%stable Variogram
%f(h) = b(1).*(1-exp(-(h./b(2))^b(3))) + b(4)
% b(1): sill (nugget removed)
% b(2): range parameter 
% b(3): smooth parameter alpha
% b(4): nugget
% Note that the nugget effect will be removed if nugget==false
%
%Matern Variogram 
% b(1): sill (nugget removed)
% b(2): range parameter
% b(3): smoothness parameter 
% b(4): nugget
%
% When nugget = false, b(4) = 0.
%
%Covariance function is given by: C(h) = C(0) - f(h) 
%where C(0) is the variance of the process and f(h) is the variogram

options = optimoptions('lsqcurvefit', 'Display', 'off');

% basic for when nugget not considered
if strcmp(fitype,'matern')
        %Matern model - no nugget
        f_base = @(b,xdata) b(1).*(1 - (1./(2.^(b(3)-1) .* gamma(b(3)))) .*...
                        ((2*sqrt(b(3)).*xdata)./b(2)).^b(3) .*...
                          besselk(b(3),(2*sqrt(b(3)).*xdata)./b(2)) );
        f = @(b,xdata) b(1).*(1 - (1./(2.^(b(3)-1) .* gamma(b(3)))) .*...
                        ((2*sqrt(b(3)).*xdata)./b(2)).^b(3) .*...
                          besselk(b(3),(2*sqrt(b(3)).*xdata)./b(2)) ) + b(4);
        fcov = @(b,xdata)(xdata>0).*( b(1)  -  b(1).*(1 - (1./(2.^(b(3)-1) .* gamma(b(3)))) .*...
                                ((2*sqrt(b(3)).*xdata)./b(2)).^b(3) .*special_besselk(b(3),(2*sqrt(b(3)).*xdata)./b(2)) ));
        lb = [0,0,0];
        ub = [inf,inf,inf];
    elseif strcmp(fitype,'stable')
        f_base = @(b, xdata) b(1) .* (1 - exp(-(xdata./b(2)).^b(3)));
        f = @(b, xdata) b(1) .* (1 - exp(-(xdata./b(2)).^b(3))) + b(4);
        fcov = @(b,xdata)(xdata>0).*( b(1)  -  b(1) .* (1 - exp(-(xdata./b(2)).^b(3))));
        lb = [0,0,0];
        ub = [2*max(v),inf,2];
    else
        error('invalid fitype %s, must be matern or stable',fitype);
end


if isempty(PrecomputeVariance)
    PrecomputeVariance = max(v);
end
% x0 = [PrecomputeVariance,1,1];
x0 = [PrecomputeVariance,min(h),1];

if ~nugget
    b=lsqcurvefit(f_base,x0,h,v,lb,ub,options);
    b = [b,0];
else
    x0 = [x0,0];
    lb = [lb,0];
    ub = [ub,0.5*PrecomputeVariance];
    b=lsqcurvefit(f,x0,h,v,lb,ub,options);
end
c_para=fcov(b,full(d_mat)); 
c_para(logical(eye(size(c_para))))=b(1) + b(4); %set diagonal to variance

end 