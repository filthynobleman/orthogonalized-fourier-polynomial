function part = voronoi_decomposition(M, sample, varargin)
%VORONOI_DECOMPOSITION(M, sample) This function computes the Voronoi
%decomposition of the mesh in the regions identified by the sample.
%
%   This function gets in input a mesh and a list of vertex indices.
%   It computes the Voronoi decomposition of the mesh into the regions
%   determined by the given sample. Each vertex of the sample identifies
%   the "center" of a region.
%   The decomposition is returned as a column vector of length n, where n
%   is the number of vertices in the mesh, and where in each row i there is
%   the index of the Voronoi region at which the i-th vertex of the mesh
%   belongs.
%   Notice that the index does not necessarly corresponds to the index of
%   the vertex in the sample, but to its position in the sample array.
%
%
%VORONOI_DECOMPOSITION(M, sample, dist_fun) This function computes the
%Voronoi decomposition of the mesh in the regions identified by the sample,
%using the given metric.
%
%   This function gets in input a mesh, a list of vertex indices and a
%   distance function.
%   It computes the Voronoi decomposition of the mesh into the regions
%   determined by the given sample. Each vertex of the sample identifies
%   the "center" of a region. The neighborhood of a sample point is
%   identified using the given distance function.
%   The decomposition is returned as a column vector of length n, where n
%   is the number of vertices in the mesh, and where in each row i there is
%   the index of the Voronoi region at which the i-th vertex of the mesh
%   belongs.
%   Notice that the index does not necessarly corresponds to the index of
%   the vertex in the sample, but to its position in the sample array.

    % Get the distance function
    dist_fun = @(Mesh, vert) (mesh.metrics.lp_dist(Mesh, vert));
    if nargin == 3
        dist_fun = varargin{1};
    end
    
    % Compute the distance matrix. In each cell (i, j) there is the
    % distance from the i-th vertex of the mesh from the j-th vertex of the
    % sample
    dists = zeros(M.n, length(sample));
    for i = 1:length(sample)
        dists(:, i) = dist_fun(M, sample(i));
    end
    
    % For each row, get the index for which the distance is minimized
    [~, part] = min(dists, [], 2);
end

