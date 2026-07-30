function [vertices, data1, data2] = prepare_data_in_double(mesh_path, map1_path, map2_path)
%This function prepare surf.gii (vertex coordinates) and shape.gii (data values) 
% file into double datatype that is required for model fitting.
%  mesh_path: path to the mesh surf.gii file
%  map1_path: path to the first map shape.gii file
%  map2_path: path to the second map

mesh_path = char(mesh_path);
map1_path = char(map1_path);
map2_path = char(map2_path);

mesh = gifti(mesh_path);
vertices = double(mesh.vertices);

map1 = gifti(map1_path);
data1 = double(map1.cdata);
map2 = gifti(map2_path);
data2 = double(map2.cdata);

end