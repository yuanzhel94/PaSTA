function [parc_idx, unique_parcs, n_parc] = get_parc(parc_arg, coord, n_points_from_range, setting, map_idx)
%[parc_idx, unique_parcs, n_parc] = get_parc(parc_arg, n_points_from_range, setting, map_idx)
%This function get the parcellation for a map based on user-specification,
%data coordinates, and global variogram fit,

if isempty(parc_arg)
    parc_idx = ones(setting.N, 1);
    unique_parcs = [];
    n_parc = 1;
else
    n_points = max(n_points_from_range,setting.min_cluster_size);
    n_parc = min(floor(setting.N / n_points),setting.max_clusters); % less than max clusters
    n_parc = max(n_parc,setting.min_clusters); % more than min clusters
    
    if isnumeric(parc_arg)
        if length(parc_arg) ~= setting.N
            error('user specified parcellation has a size %d differ from data %d', length(parc_arg), setting.N);
        end
        [unique_parcs, ~, parc_idx] = unique(parc_arg, 'sorted');
        n_parc_in = length(unique_parcs);
        if n_parc_in > n_parc
            warning(['data No.%d: specified number of parcs %d is larger than data-derived max number of parcs %d, ' ...
                'carefully tradeoff the ability for detecting nonstationary and the parcel coverage for robust estimation'], map_idx, n_parc_in, n_parc);
        end
        n_parc = n_parc_in;

    elseif strcmp(parc_arg,'auto')
        if n_parc == 1
            parc_idx = ones(setting.N, 1);
            unique_parcs = 1;
        else
            parc_idx = kmeans(coord, n_parc);
            [unique_parcs, ~, parc_idx] = unique(parc_idx, 'sorted');
        end
        
    else
        error('invalid parcellation input type for map %d', map_idx);
    end
end


end