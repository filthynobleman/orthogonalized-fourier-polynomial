function M = compute_laplace_basis(M, k)
%M = compute_laplace_basis(M, k) Computes the first k eigenfunctions of the
%Laplace-Beltrami basis for the given mesh. Notice that the first constant
%eigenfunction is excluded from the basis.
%
%   This is a full list of the newly computed fields of the structure:
%       - A:            The lumped mass matrix. If present, it is not
%                       recomputed.
%       - S:            The stiffness matrix. If present, it is not
%                       recomputed.
%       - evecs:        The first k + 1 eigenfuntions of the
%                       Laplace-Beltrami basis, including the constant one.
%       - evals:        The eigenvalues associated to 'evecs'.
%       - basis:        The first k eigenfunctions of the Laplace-Beltrami 
%                       basis, excluding the constant one.
%       - support:      The eigenvalues associated to 'basis'.
%       - phi0:         The value assumed by the constant eigenfunction in
%                       the whole mesh.
%       - basis_len:    The length of the basis. Namely, k.

    if ~isfield(M, 'A') || ~isfield(M, 'S')
        [M.S, M.A, M.Al] = mesh.proc.laplacian(M);
    end
    if isfield(M, 'basis') && size(M.basis, 2) >= k
        return;
    end
    if k >= M.n
        warning(strcat("The number of eigenfunctions (%d) exceeds ", ...
                       "rank of the laplacian (%d)."), k, M.n - 1);
        k = M.n - 1;
    end
    [M.evecs, M.evals] = eigs(M.S, M.A, k+1, -1e-5);
    M.evals = diag(M.evals);
    M.support = M.evals(2:k+1);
    M.basis = M.evecs(:, 2:k+1);
    M.phi0 = M.evecs(1, 1);
    M.basis_len = k;
end

