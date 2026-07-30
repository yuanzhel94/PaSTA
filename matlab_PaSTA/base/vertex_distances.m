function [ed, gd] = vertex_distances(v,f)
%Compute the Euclidean space (ed) and geodesic distances (gd) between all
%pairs of vertices. Geodesic distances are computed by defining adjacency
%graph between vertices, with edge weight being Euclidean distance between
%connected vertices and 0 between disconnected vertices.
%
%v:
% N-by-3 matrix of vertex coordinates, each row is the coordinate of a vertex.
%f:
% M-by-3 faces definition, each row is the index of three vertices defining
% a face.

ed = squareform(pdist(v));
adj = sparse(size(v,1), size(v,1));

for i = 1:size(f, 1)
    face = f(i, :);
    for j = 1:3
        k = mod(j, 3) + 1;
        d = ed(face(j),face(k));
        adj(face(j),face(k)) = d;
        adj(face(k),face(j)) = d;
    end
end

G = graph(adj);
gd = distances(G);

end