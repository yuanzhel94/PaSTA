function [weight_sum, weighted_diff2_sum_x] = blockij_variogram_accumulation_single_map(d, xi, xj, spec, is_diagonal_block)
%[weight_sum, weighted_diff2_sum_x] = blockij_variogram_accumulation_single_map(d, xi, xj, spec, is_diagonal_block)
%This function computes the accumulated weights and weighted difference
%squares for block ij for a single map. These are later used to compute a global variogram.

xi = xi(:);
xj = xj(:);
h = spec.h(:);

weight_sum = zeros(spec.M, 1);
weighted_diff2_sum_x = zeros(spec.M, 1);

% pairs fail within variogram evaluation range; triu for diagonal blocks
if is_diagonal_block
    valid = triu((d > 0) & (d <= spec.variogram_dmax), 1);
else
    valid = (d > 0) & (d <= spec.variogram_dmax);
end

if ~any(valid, 'all') % in case no pair in variogram range
    return
end

d = d(valid);
diff2_x = (xi - xj.') .^ 2;
diff2_x = diff2_x(valid);
d = d(:);
diff2_x = diff2_x(:);

nearest_bin = round((d - spec.variogram_dmin) / spec.lag_sep) + 1;
for offset = spec.bin_idx_offset
    bin_idx = nearest_bin + offset;
    valid_bins = (bin_idx >= 1) & (bin_idx <= spec.M);
    if ~any(valid_bins) % allocated to bins
        continue;
    end
    valid_bin_idx = bin_idx(valid_bins);
    d_valid_bins = d(valid_bins);
    valid_bin_idx = valid_bin_idx(:);
    d_valid_bins = d_valid_bins(:);
    d2_from_bins = (d_valid_bins - h(valid_bin_idx)) .^ 2;
    inside_kernel = d2_from_bins <= spec.L2;
    if ~any(inside_kernel) % within the cutoff kernel size
        continue;
    end
    valid_bin_idx = valid_bin_idx(inside_kernel);
    d2_from_bins = d2_from_bins(inside_kernel);
    final_idx = find(valid_bins);
    final_idx = final_idx(inside_kernel); % idx of pairs that contribute to final weights
    weight = exp(spec.weight_constant .* d2_from_bins);
    weight_sum = weight_sum + accumarray(valid_bin_idx, weight, [spec.M, 1], @sum, 0); % accumulate and pad into shape (M, 1)
    weighted_diff2_sum_x = weighted_diff2_sum_x + accumarray(valid_bin_idx, weight .* diff2_x(final_idx), [spec.M, 1], @sum, 0);
end

end