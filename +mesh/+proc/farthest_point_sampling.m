function fps = farthest_point_sampling(M, n, varargin)
%FARTHEST_POINT_SAMPLING(M, n) This function generates a farthest point
%sampling of n vertices from the given mesh.
%
%   This function gets in input a mesh and an integer value.
%   It returns a column array of length n containing the indices of the
%   vertices of the mesh which are part of the farthest point sampling. The
%   farthest point sampling is computed using the Euclidean L2 distance as
%   metric.
%
%
%FARTHEST_POINT_SAMPLING(M, n, dist_fun) This function generates a farthest
%point sampling of n vertices from the given mesh.
%
%   This function gets in input a mesh and, integer value and a distance
%   function.
%   It returns a column array of length n containing the indices of the
%   vertices of the mesh which are part of the farthest point sampling. The
%   farthest point sampling is computed using the given distance function
%   as metric.

    % Get the distance function, if given
    dist_fun = @(Mesh, vert) (mesh.metrics.lp_dist(Mesh, vert));
    if nargin == 3
        dist_fun = varargin{1};
    end
    
    % Initialize the FPS with a random vertex
    v = 1;%uint64(rand(1) * M.n);
    fps = zeros(n, 1);
    fps(1) = v;
    for k = 1:(n - 1)
        dists = set_dist(M, fps(1:k), dist_fun);
        [~, idxs] = max(dists);
        fps(k + 1) = idxs(1);
    end
end

function dists = set_dist(M, set, dist_fun)
    dists = dist_fun(M, set(1));
    for i = 2:max(size(set))
        dists = min(dists, dist_fun(M, set(i)));
    end
end