function [nef, run_status] = covs2nef(cov1,cov2)
% efficient implementation of effective sample size computation from
% covariance matrices
cov1(~isfinite(cov1)) = 0;
cov2(~isfinite(cov2)) = 0;
c1 = cov1 - mean(cov1,1) - mean(cov1,2) + mean(cov1(:));
c2 = cov2 - mean(cov2,1) - mean(cov2,2) + mean(cov2(:));
num = trace(c1 * c2);
den = trace(c1) * trace(c2);
nef = real(1 / (num / den) + 1);
run_status = nef > 2;

% equivalent to the computation below but more computationally efficient:
% nef=real(1/(trace(B*fc_para1*B*fc_para2)/(trace(B*fc_para1)*trace(B*fc_para2)))+1);

end