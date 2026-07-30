function [parc, n_parc, unique_parcs, fc_para] = parc_data(parc, c_para, b, D, coord, max_clusters, min_clusters, min_cluster_size, map_idx)
% parcellating data if needed
% fc_para will be returned same as c_para if only 1 parc, otehrwise a zero
% matrix with the same shape of c_para

if isempty(parc)
    fc_para = c_para;
    n_parc = 1;
    parc = [];
    unique_parcs = [];
else
    range = b(2) * 2.996 ^ (1/b(3));
    n_points = sum(D(:) < range) / size(D,1) - 1;
    n_points = max(n_points,min_cluster_size);
    n_clusters = min(floor(size(D,1) / n_points),max_clusters); % less than max clusters
    n_clusters = max(n_clusters,min_clusters); % more than min clusters
    if strcmp(parc,'auto')
        parc = kmeans(coord,n_clusters);
        unique_parcs = unique(parc);
        n_parc = length(unique_parcs);
    else
        unique_parcs = unique(parc);
        n_parc = length(unique_parcs);
        if n_parc > n_clusters
            warning(['data No.%d: specified number of parcs %d is larger than data-derived max number of parcs %d, ' ...
                'carefully tradeoff the ability for detecting nonstationary and the parcel coverage for robust estimation'], map_idx, n_parc, n_clusters);
        end
    end
    if n_parc == 1
        fc_para = c_para;
    else
        fc_para = zeros(size(c_para));
    end
end

end