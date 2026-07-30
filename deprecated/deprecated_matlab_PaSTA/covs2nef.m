function [nef, run_status] = covs2nef(cov1,cov2)
% efficient implementation of effective sample size computation from
% covariance matrices
cov1(~isfinite(cov1)) = 0;
cov2(~isfinite(cov2)) = 0;

N = size(cov1, 1);
r1 = sum(cov1, 2);
r2 = sum(cov2, 2);
s1 = sum(r1);
s2 = sum(r2);
tr1 = sum(diag(cov1));
tr2 = sum(diag(cov2));
A1 = tr1 - s1 / N;
A2 = tr2 - s2 / N;
tr12 = sum(cov1 .* cov2, 'all');
num = tr12 - 2 * (r1' * r2) / N + (s1 * s2) / N^2;
den = A1 * A2;
nef = real(den / num + 1);

run_status = isfinite(nef) && isfinite(num) && num > 0 && nef > 2;

% c1 = cov1 - mean(cov1,1) - mean(cov1,2) + mean(cov1(:));
% c2 = cov2 - mean(cov2,1) - mean(cov2,2) + mean(cov2(:));
% num = trace(c1 * c2);
% den = trace(c1) * trace(c2);
% nef = real(1 / (num / den) + 1);
% run_status = nef > 2;

% equivalent to the computation below but more computationally efficient:
% nef=real(1/(trace(B*fc_para1*B*fc_para2)/(trace(B*fc_para1)*trace(B*fc_para2)))+1);

end