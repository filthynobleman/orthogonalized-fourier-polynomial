function idx = unique_indices(k, order)
%%WARNING: assuming order 2 for now
if nargin < 2
    order = 2;
end
idx = [1:k];%zeros((k^2 + k)/2, 1);
for i = 1:k
    idx = [idx (i*k + i):((i + 1)*k)];
end
idx = [1, (idx + 1)];
end

