function K = heat_kernel(M, t, k)
%K = HEAT_KERNEL(M, t, k) Computes the approximation of the heat kernel of
%a mesh with k eigenfunctions of the Laplacian.
%
%K = HEAT_KERNEL(M, t) Computes the approximation of the heat kernel of the
%mesh with all the eigenfunctions already computed.
%
%   It computes the approximation of the heat kernel of the mesh M, using
%   the first k eigenfunctions of the laplacian of M. The computed heat
%   kernel is parametric with respect to time instant t.
%   The heat kernel is represented as an n-by-n matrix, where n is the
%   number of vertices of M and each cell (i, j) contains the heat
%   diffusion from the vertex i from vertex j.

    % If the laplace basis has not yet been computed, compute it
    if ~isfield(M, 'basis')
        if nargin < 3
            error(strcat("If the given mesh has not a precomputed ", ...
                         "basis, then a number of eigenfunctions for ", ...
                         "the approximation must be given."));
        end
        M = mesh.proc.compute_laplace_basis(M, k);
    % If the laplace basis is not large enough, the recompute it
    elseif nargin == 3 && size(M.basis, 2) < k
        M = mesh.proc.compute_laplace_basis(M, k);
    % Otherwise, a basis is provided, but not approximation limit is given.
    % In this case, use all of the available basis
    else
        k = size(M.basis, 2);
    end
    % Compute the heat kernel
    K = M.basis(:, 1:k) * diag(exp(- M.support * t)) * M.basis(:, 1:k)';
end

