function energies = dirichlet(M, f)
%energies = DIRICHLET(M, f) Computes the Dirichlet energy of the given
%function. If f is a matrix, it is treated as a set of functions, one for
%each column.
    
    % Check if M has a mass and a stiffness matrix
    if ~isfield(M, 'S') || ~isfield(M, 'A')
        [M.S, ~, M.A] = mesh.proc.laplacian(M);
    end
    % Compute the dirichlet energies
    energies = diag(f' * (M.S * f))';
end

