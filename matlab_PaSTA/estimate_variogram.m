function [v,h] = estimate_variogram(D, data, M, qd)
% Estimate the empirical variogram from distance matrix between vertices,
% and data value at each vertex. Estiamtion performed in M bins, ranging
% from min_distance to qd * max_distance, where max_distance is the max distance in
% the distance matrix.
%
%D:
% Distance matrix between all vertices
%data:
% Data value at each vertex.
%M:
% number of bins to estimate variogram.
%qd:
% maximum distance to evaluate variogram
%
% Returns:
% v: semivariance of the varigoram (M,1)
% h: lag distances of the variogram (1, M)
Dmax = qd * max(D(:));
Dmin = min(D(D~=0));
triu_D = triu(D,1);
[row,col,dval]=find((triu_D <= Dmax) .* triu_D);
h=linspace(Dmin,Dmax,M);
delta=(Dmax-Dmin)/(M-1)*0.5;
sigma=6*delta;
v=zeros(M,1);
for i=1:M
    w=exp( -1*((2.68*abs(h(i)-dval)).^2)/(2*sigma^2) ); 
    v(i)=0.5 * sum(w.* ( data(row) - data(col) ).^2); 
    v(i)=v(i)/sum(w); 
end

end