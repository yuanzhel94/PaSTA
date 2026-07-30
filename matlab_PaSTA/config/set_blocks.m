function block_setting = set_blocks(N, block_size, n_workers)
%block_setting = set_blocks(N, block_size, n_workers)
%This function specifies block setting used for memory efficiency
%

N = double(N);
n_workers = double(n_workers);

if isempty(block_size)
    block_size = N;
else
    block_size = min(double(block_size), N);
end

block_setting = struct();
block_setting.N = N;
block_setting.block_size = block_size;
block_setting.n_workers = n_workers;

block_setting.n_blocks = ceil(N / block_size);
block_setting.block_start = (0:(block_setting.n_blocks - 1))' .* block_size + 1;
block_setting.block_end = min(block_setting.block_start + block_size - 1, N);

block_setting.n_unique_block_pairs = block_setting.n_blocks * (block_setting.n_blocks + 1) / 2;
block_setting.block_i = zeros(block_setting.n_unique_block_pairs, 1);
block_setting.block_j = zeros(block_setting.n_unique_block_pairs, 1);

pointer_idx = 1;
for block_i = 1:block_setting.n_blocks
    n_block_j = block_setting.n_blocks - block_i + 1;
    pair_idx = pointer_idx:(pointer_idx + n_block_j - 1);
    block_setting.block_i(pair_idx) = block_i;
    block_setting.block_j(pair_idx) = (block_i:block_setting.n_blocks)';
    pointer_idx = pointer_idx + n_block_j;
end

block_setting.is_diagonal_block = block_setting.block_i == block_setting.block_j;

block_setting.parallel = n_workers > 1 && block_setting.n_unique_block_pairs > 1;

end