function distance_memory = prepare_distance_in_memory(distance_data, N)
%distance_memory = prepare_distance_in_memory(distance_data, N)
%This function validates and prepares an in-memory distance input.
%
%distance_data can be:
%   1. A full square distance matrix of shape (N, N).
%   2. A vector containing the strict upper-triangular elements in MATLAB
%      column-major order:
%      D(1,2), D(1,3), D(2,3), D(1,4), D(2,4), D(3,4), ...
%
%The input representation and numeric precision are retained. Conversion
%to double is delayed until individual distance blocks are extracted.

if ~isnumeric(distance_data) || ~isreal(distance_data) || issparse(distance_data)
    error('distance_data must be a real, nonsparse numeric array.');
end

if ~isnumeric(N) || ~isscalar(N) || ~isreal(N) || ...
        ~isfinite(N) || N < 1 || fix(N) ~= N
    error('N must be a positive integer scalar.');
end

N = double(N);
n_upper = N * (N - 1) / 2;

if ismatrix(distance_data) && ...
        size(distance_data, 1) == N && ...
        size(distance_data, 2) == N

    storage = 'matrix';

elseif isvector(distance_data) && numel(distance_data) == n_upper

    storage = 'upper';

else
    error(['distance_data must be either an N-by-N distance matrix or ' ...
        'a strict upper-triangular vector containing ' ...
        'N * (N - 1) / 2 = %.0f elements.'], n_upper);
end

distance_memory = struct();
distance_memory.data = distance_data;
distance_memory.storage = storage;
distance_memory.N = N;
distance_memory.n_upper = n_upper;
distance_memory.data_class = class(distance_data);

end