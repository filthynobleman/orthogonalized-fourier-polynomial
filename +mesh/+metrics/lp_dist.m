function dists = lp_dist(varargin)
%dists = LP_DIST(v1, v2)
%
%   Computes the L2 distance between vertex v1 and vertex v2. If v1 and v2
%   are matrices, they must be of the same length and the function computes
%   the distance row-wise and returns a vector containing in the i-th 
%   position the distance between v1(i) and v2(i).
%
%
%dists = LP_DIST(M, v)
%
%   Computes the Euclidean L2 distance from the given vertex to all the
%   vertices of the mesh. The vertex v can be given as an integer index of
%   a vertex in the mesh M or as a triplet of real coordinates.
%
%
%dists = LP_DIST(--, p)
%
%   Computes the Euclidean Lp distance, rather than the L2 one, according
%   to the possible parameters of the function.

    % Check if p is given or if it is two
    if nargin == 3
        p = varargin{3};
    else
        p = 2;
    end
    
    % Check if the first argument is a mesh or a vertex/matrix
    V1 = varargin{1};
    V2 = varargin{2};
    if max(size(V1)) == 1
        V1 = V1.VERT;
        if max(size(V2)) == 1
            V2 = V1(V2, :);
        end
    end

    % Compute the distance
    dists = sum(abs(V1 - V2).^p, 2).^(1/p);
end

