function [vx, vy] = parcel_variogram_with_blocks_loaded(x, y, distance_matfile, distance_variable, parc_idx_x, parc_idx_y, parcel_specs_x, parcel_specs_y, distance_blocks)
%[vx, vy] = parcel_variogram_with_blocks_loaded(x, y, distance_matfile, distance_variable, parc_idx_x, parc_idx_y, parcel_specs_x, parcel_specs_y, distance_blocks)
%This function computes parcel-level variograms for two maps using loaded distance
%matrix in blockwise manner

parc_idx_x = parc_idx_x(:);
parc_idx_y = parc_idx_y(:);

n_parc_x = numel(parcel_specs_x);
n_parc_y = numel(parcel_specs_y);

Mx = cellfun(@(spec) spec.M, parcel_specs_x);
My = cellfun(@(spec) spec.M, parcel_specs_y);

Mx = Mx(:);
My = My(:);

end_idx_x = cumsum(Mx);
start_idx_x = end_idx_x - Mx + 1;

end_idx_y = cumsum(My);
start_idx_y = end_idx_y - My + 1;

total_M_x = sum(Mx);
total_M_y = sum(My);

weight_sum_x = zeros(total_M_x, 1);
weighted_diff2_sum_x = zeros(total_M_x, 1);

weight_sum_y = zeros(total_M_y, 1);
weighted_diff2_sum_y = zeros(total_M_y, 1);

if distance_blocks.parallel
    parfor (pair_idx = 1:distance_blocks.n_unique_block_pairs, distance_blocks.n_workers)

        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        is_diagonal_block = distance_blocks.is_diagonal_block(pair_idx);

        [dij, idx_i, idx_j] = load_valid_distance_blockij(distance_matfile, ...
            distance_variable, ...
            block_i, ...
            block_j, ...
            distance_blocks);

        weight_sum_x_ij = zeros(total_M_x, 1);
        weighted_diff2_sum_x_ij = zeros(total_M_x, 1);

        weight_sum_y_ij = zeros(total_M_y, 1);
        weighted_diff2_sum_y_ij = zeros(total_M_y, 1);

        if ~isempty(dij)
            [weight_sum_x_cell, weighted_diff2_sum_x_cell] = ...
                blockij_parcel_variogram_accumulation_single_map(dij, x(idx_i), x(idx_j), parc_idx_x(idx_i), parc_idx_x(idx_j), parcel_specs_x, is_diagonal_block);

            [weight_sum_y_cell, weighted_diff2_sum_y_cell] = ...
                blockij_parcel_variogram_accumulation_single_map(dij, y(idx_i), y(idx_j), parc_idx_y(idx_i), parc_idx_y(idx_j), parcel_specs_y, is_diagonal_block);

            for parc = 1:n_parc_x
                selected = start_idx_x(parc):end_idx_x(parc);
                weight_sum_x_ij(selected) = weight_sum_x_cell{parc};
                weighted_diff2_sum_x_ij(selected) = weighted_diff2_sum_x_cell{parc};
            end

            for parc = 1:n_parc_y
                selected = start_idx_y(parc):end_idx_y(parc);
                weight_sum_y_ij(selected) = weight_sum_y_cell{parc};
                weighted_diff2_sum_y_ij(selected) = weighted_diff2_sum_y_cell{parc};
            end
        end

        weight_sum_x = weight_sum_x + weight_sum_x_ij;
        weighted_diff2_sum_x = weighted_diff2_sum_x + weighted_diff2_sum_x_ij;

        weight_sum_y = weight_sum_y + weight_sum_y_ij;
        weighted_diff2_sum_y = weighted_diff2_sum_y + weighted_diff2_sum_y_ij;
    end

else
    for pair_idx = 1:distance_blocks.n_unique_block_pairs

        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);
        is_diagonal_block = distance_blocks.is_diagonal_block(pair_idx);

        [dij, idx_i, idx_j] = load_valid_distance_blockij(distance_matfile, ...
            distance_variable, ...
            block_i, ...
            block_j, ...
            distance_blocks);

        weight_sum_x_ij = zeros(total_M_x, 1);
        weighted_diff2_sum_x_ij = zeros(total_M_x, 1);

        weight_sum_y_ij = zeros(total_M_y, 1);
        weighted_diff2_sum_y_ij = zeros(total_M_y, 1);

        if ~isempty(dij)
            [weight_sum_x_cell, weighted_diff2_sum_x_cell] = ...
                blockij_parcel_variogram_accumulation_single_map(dij, x(idx_i), x(idx_j), parc_idx_x(idx_i), parc_idx_x(idx_j), parcel_specs_x, is_diagonal_block);

            [weight_sum_y_cell, weighted_diff2_sum_y_cell] = ...
                blockij_parcel_variogram_accumulation_single_map(dij, y(idx_i), y(idx_j), parc_idx_y(idx_i), parc_idx_y(idx_j), parcel_specs_y, is_diagonal_block);
            

            for parc = 1:n_parc_x
                selected = start_idx_x(parc):end_idx_x(parc);
                weight_sum_x_ij(selected) = weight_sum_x_cell{parc};
                weighted_diff2_sum_x_ij(selected) = weighted_diff2_sum_x_cell{parc};
            end

            for parc = 1:n_parc_y
                selected = start_idx_y(parc):end_idx_y(parc);
                weight_sum_y_ij(selected) = weight_sum_y_cell{parc};
                weighted_diff2_sum_y_ij(selected) = weighted_diff2_sum_y_cell{parc};
            end
        end

        weight_sum_x = weight_sum_x + weight_sum_x_ij;
        weighted_diff2_sum_x = weighted_diff2_sum_x + weighted_diff2_sum_x_ij;

        weight_sum_y = weight_sum_y + weight_sum_y_ij;
        weighted_diff2_sum_y = weighted_diff2_sum_y + weighted_diff2_sum_y_ij;
    end
end

vx = cell(n_parc_x, 1);
vy = cell(n_parc_y, 1);

for parc = 1:n_parc_x
    selected = start_idx_x(parc):end_idx_x(parc);
    weight_sum_parc = weight_sum_x(selected);
    weighted_diff2_sum_parc = weighted_diff2_sum_x(selected);
    vx{parc} = nan(Mx(parc), 1);
    valid_h = weight_sum_parc > 0;
    vx{parc}(valid_h) = 0.5 * weighted_diff2_sum_parc(valid_h) ./ weight_sum_parc(valid_h);
end

for parc = 1:n_parc_y
    selected = start_idx_y(parc):end_idx_y(parc);
    weight_sum_parc = weight_sum_y(selected);
    weighted_diff2_sum_parc = weighted_diff2_sum_y(selected);
    vy{parc} = nan(My(parc), 1);
    valid_h = weight_sum_parc > 0;
    vy{parc}(valid_h) = 0.5 * weighted_diff2_sum_parc(valid_h) ./ weight_sum_parc(valid_h);
end

end