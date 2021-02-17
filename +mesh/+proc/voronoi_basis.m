function varargout = voronoi_basis(varargin)
%V = VORONOI_BASIS(M, k)
%
%   Computes the k-th dimensional Voronoi basis for the mesh M using the
%   Euclidean L2 distance function.
%
%
%V = VORONOI_BASIS(M, k, voronoi_dist)
%
%   Computes the k-dimensional Voronoi basis for the mesh M using
%   voronoi_dist as distance function and the Euclidean L2 distance as
%   sampling distance function.
%
%
%V = VORONOI_BASIS(M, k, voronoi_dist, fps_dist)
%
%   Computes the k-dimensional Voronoi basis for the mesh M using
%   voronoi_dist as distance function and sampling the mesh using fps_dist
%   as distance function.
%
%
%V = VORONOI_BASIS(M, sample)
%
%   Computes the Voronoi basis for the mesh M using the given sample ans
%   the Euclidean L2 distance function. The dimension of the basis is the
%   same as the length of the sample.
%
%
%V = VORONOI_BASIS(M, sample, voronoi_dist)
%
%   Computes the Voronoi basis for the mesh M using the given distance
%   function. If the sample is not given, the same distance function is
%   also used for sampling the mesh.
%
%
%[V, orthoV] = VORONOI_BASIS(--)
%
%   Computes the Voronoi basis and its normalized version, which is
%   orthonormal.

    % Get the mesh
    M = varargin{1};
    % If the sample is given, use it and set k to its length
    if length(varargin{2}) > 1
        sample = varargin{2};
        k = length(sample);
    % Otherwise, set k to the second parameter and perform the sampling
    else
        k = varargin{2};
        % If the distance function is the default one, sample with default
        if nargin < 4
            sample = mesh.proc.farthest_point_sampling(M, k);
        % Otherwise, use the given one
        else
            fps_dist = varargin{4};
            sample = mesh.proc.farthest_point_sampling(M, k, fps_dist);
        end
    end
    
    % Compute the Voronoi decomposition with the proper distance function
    if nargin < 3
        voronoi = mesh.proc.voronoi_decomposition(M, sample);
    else
        voronoi_dist = varargin{3};
        voronoi = mesh.proc.voronoi_decomposition(M, sample, voronoi_dist);
    end
    
    % Compute the Voronoi basis
    V = zeros(M.n, k);
    for i = 1:k
        V(:, i) = double(voronoi == i);
    end
    
    % Returns the first output value
    varargout{1} = V;
    % If the orthonormal version is required, compute it
    if nargout == 2
        orthoV = V;
        for i = 1:k
            orthoV(:, i) = V(:, i) / norm(V(:, i));
        end
        varargout{2} = orthoV;
    end
    
end

