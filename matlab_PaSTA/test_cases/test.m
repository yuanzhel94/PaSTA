map0_path = fullfile('fsaverage10k', '0.shape.gii');
map1_path = fullfile('fsaverage10k', '1.shape.gii');
mesh_path = fullfile('fsaverage10k', 'mesh.surf.gii');
[coords, x, y] = prepare_data_in_double(mesh_path, map0_path, map1_path);

D = squareform(pdist(coords));
distance_file_path = fullfile('fsaverage10k', 'D_full.mat');
save(distance_file_path, 'D', '-v7.3');

N = size(D,1);
D_triu = single(D(triu(true(N), 1)));

tic; pasta_fit = pasta(x, y, coords); toc;
tic; pasta_fit = pasta(x, y, coords, 'xparc', 'auto', 'yparc', 'auto', 'random_state', 0); toc;

tic; pasta_fit = pasta(x, y, coords, 'D', D); toc;
tic; pasta_fit = pasta(x, y, coords, 'D', D, 'xparc', 'auto', 'yparc', 'auto', 'random_state', 0); toc;

tic; pasta_fit = pasta(x, y, coords, 'D', D_triu); toc;
tic; pasta_fit = pasta(x, y, coords, 'D', D_triu, 'xparc', 'auto', 'yparc', 'auto', 'random_state', 0); toc;

tic; pasta_fit = pasta(x, y, coords, 'distance_file', distance_file_path, 'distance_variable', 'D'); toc;
tic; pasta_fit = pasta(x, y, coords, 'distance_file', distance_file_path, 'distance_variable', 'D', 'xparc', 'auto', 'yparc', 'auto', 'random_state', 0); toc;
