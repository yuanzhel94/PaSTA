function [fc_para, bs] = parcellated_cov(x, parc, unique_parcs, D, dim, b, fc_para, qd, M_scale, nugget, crossblock)
%This functon computes nonstationary covariance matrix of data x based on
%its parcellation provided (parc)
% Inputs:
%  x: spatial map data, vector of length N
%  parc: membership of parcel for each data in x, int vector of length N
%  unique_parcs: unique ints in parc each corresponding to a parcel
%  D: distance matrix (N,N) between pairwise data in x
%  dim: dimension of the spatial coordinates
%  b: statioanry stable variogram model parameter previously fitted using all data
%  fc_para: initialized covariance matrix
%  qd: qd determines the distance used when estimating variogram
%  M_scale: scaling factor of the number of lag distances evaluated
%  nugget: bool indicator of whether using nugget, should use true
%  crossblock: ways to infer cross-block covariance, always use 'conv'
%       (i.e., process convolution) despite 'independent' (assuming cross-parcel
%       independence) is also implemented.
%
% Returns:
%   fc_para: fitted nonstationary covariance matrix
%   bs: fitted stable variogram model parameter for each parcel
exponent = b(3);
n_parc = length(unique_parcs);
bs = zeros(n_parc, 4);
% fit stable variogram and covariance within each parcel using
% predetermined shape parameter
for i=1:n_parc
    parci = unique_parcs(i);
    v_select = parc==parci;
    x_select = x(v_select);var_x_select = var(x_select); x_select = zscore(x_select);
    D_select = D(v_select,v_select);
    MM = round(sqrt(length(x_select))) * M_scale;
    [va_select,h1_select] = estimate_variogram(D_select, x_select, MM, qd);
    [pc_para,pb]=fit_stable_fix_exponent(h1_select',va_select,D_select,[],exponent,nugget); pb(1) = pb(1) * var_x_select;pb(4) = pb(4) * var_x_select;
    fc_para(v_select,v_select) = pc_para .* var_x_select;
    bs(i,:) = pb;
end

% infer cross-parcel covariance, using either independence assumption
% between parcels (i.e., independent, not realistic) or process convolution
% (i.e., 'conv')
switch crossblock
    case 'independent'
        % do nothing with cross block elements, assuming between parcels
        % are independent, which is not realistic particularly for data
        % points at parcel edges
    case 'conv'
        for i=1:n_parc-1
            parci = unique_parcs(i);
            v_select1 = parc==parci;
            phii=bs(i,2);
            for j=i:n_parc
                parcj = unique_parcs(j);
                v_select2 = parc==parcj;
                D_select = D(v_select1,v_select2);
                phij=bs(j,2);
                Sigma = (phii .^ 2 + phij .^ 2)/2; 
                Qij = D_select .^ 2 ./ Sigma; 
                fc_para(v_select1,v_select2) = (phii * phij / Sigma)^(dim/2) * sqrt(bs(i,1) * bs(j,1)) * exp(- (sqrt(Qij) .^ (exponent)));
                fc_para(v_select2,v_select1) = fc_para(v_select1,v_select2)';
            end
        end
    otherwise
        error('invalid crossblock argument');
end

end