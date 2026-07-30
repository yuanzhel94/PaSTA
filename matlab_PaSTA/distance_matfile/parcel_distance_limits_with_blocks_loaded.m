function [dmin_x, dmax_x, dmin_y, dmax_y] = parcel_distance_limits_with_blocks_loaded(distance_matfile, distance_variable, parc_idx_x, n_parc_x, parc_idx_y, n_parc_y, distance_blocks)
%[dmin_x, dmax_x, dmin_y, dmax_y] = parcel_distance_limits_with_blocks_loaded(distance_matfile, distance_variable, parc_idx_x, n_parc_x, parc_idx_y, n_parc_y, distance_blocks)
%This function computes the min and max distance for each parcel using
%blockwise manner

n_pairs = distance_blocks.n_unique_block_pairs;
dmin_x_per_pair = inf(n_parc_x, n_pairs);
dmax_x_per_pair = zeros(n_parc_x, n_pairs);
dmin_y_per_pair = inf(n_parc_y, n_pairs);
dmax_y_per_pair = zeros(n_parc_y, n_pairs);

same_parcellation = n_parc_x == n_parc_y && isequal(parc_idx_x, parc_idx_y);

if distance_blocks.parallel
    parfor (pair_idx = 1:n_pairs, distance_blocks.n_workers)
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);

        [dij, idx_i, idx_j] = load_valid_distance_blockij(distance_matfile, ...
            distance_variable, ...
            block_i, ...
            block_j, ...
            distance_blocks);

        dmin_x_ij = inf(n_parc_x, 1);
        dmax_x_ij = zeros(n_parc_x, 1);
        dmin_y_ij = inf(n_parc_y, 1);
        dmax_y_ij = zeros(n_parc_y, 1);

        if ~isempty(dij)
            [dmin_x_ij, dmax_x_ij] = blockij_parcel_distance_limits(dij, parc_idx_x(idx_i), parc_idx_x(idx_j), n_parc_x);

            if same_parcellation
                dmin_y_ij = dmin_x_ij;
                dmax_y_ij = dmax_x_ij;
            else
                [dmin_y_ij, dmax_y_ij] = blockij_parcel_distance_limits(dij, parc_idx_y(idx_i), parc_idx_y(idx_j), n_parc_y);
            end
        end

        dmin_x_per_pair(:, pair_idx) = dmin_x_ij;
        dmax_x_per_pair(:, pair_idx) = dmax_x_ij;
        dmin_y_per_pair(:, pair_idx) = dmin_y_ij;
        dmax_y_per_pair(:, pair_idx) = dmax_y_ij;
    end

else
    for pair_idx = 1:n_pairs
        block_i = distance_blocks.block_i(pair_idx);
        block_j = distance_blocks.block_j(pair_idx);

        [dij, idx_i, idx_j] = load_valid_distance_blockij(distance_matfile, ...
            distance_variable, ...
            block_i, ...
            block_j, ...
            distance_blocks);

        dmin_x_ij = inf(n_parc_x, 1);
        dmax_x_ij = zeros(n_parc_x, 1);
        dmin_y_ij = inf(n_parc_y, 1);
        dmax_y_ij = zeros(n_parc_y, 1);

        if ~isempty(dij)
            [dmin_x_ij, dmax_x_ij] = blockij_parcel_distance_limits(dij, parc_idx_x(idx_i), parc_idx_x(idx_j), n_parc_x);
            if same_parcellation
                dmin_y_ij = dmin_x_ij;
                dmax_y_ij = dmax_x_ij;
            else
                [dmin_y_ij, dmax_y_ij] = blockij_parcel_distance_limits(dij, parc_idx_y(idx_i), parc_idx_y(idx_j), n_parc_y);
            end
        end

        dmin_x_per_pair(:, pair_idx) = dmin_x_ij;
        dmax_x_per_pair(:, pair_idx) = dmax_x_ij;
        dmin_y_per_pair(:, pair_idx) = dmin_y_ij;
        dmax_y_per_pair(:, pair_idx) = dmax_y_ij;
    end
end

dmin_x = min(dmin_x_per_pair, [], 2);
dmax_x = max(dmax_x_per_pair, [], 2);
dmin_y = min(dmin_y_per_pair, [], 2);
dmax_y = max(dmax_y_per_pair, [], 2);

end