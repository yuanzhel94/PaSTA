function [nef, run_status] = stats2nef(rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum)
%[nef, run_status] = stats2nef(rowsum_x, rowsum_y, trace_x, trace_y, cov_prodsum)
%This function computes the effective sample  size from blockwise
%accumulation statistics
N = length(rowsum_x);
sum_x = sum(rowsum_x);
sum_y = sum(rowsum_y);
centered_trace_x = trace_x - sum_x / N;
centered_trace_y = trace_y - sum_y / N;
cov_prod = cov_prodsum - 2 * dot(rowsum_x, rowsum_y) / N + (sum_x * sum_y) / N ^ 2;
nef = real(centered_trace_x * centered_trace_y / cov_prod + 1);
run_status = isfinite(nef) && isfinite(cov_prod) && cov_prod > 0 && nef > 2;
end