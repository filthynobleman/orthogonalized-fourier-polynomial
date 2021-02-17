function M = compute_extbasis(M, k, order, normalize)
%M = COMPUTE_EXTBASIS(M) Compute the extended basis of mesh M, containing
%the polynomials of order 2 of all the Laplace eigenfunctions of M.
%
%M = COMPUTE_EXTBASIS(M, k) Compute the extended basis of mesh M,
%containing the polynomials of order 2 of the first k Laplace 
%eigenfunctions of M.
%
%M = COMPUTE_EXTBASIS(M, k, order) Compute the extended basis of mesh M,
%containing the polynomials of the given order of the first k Laplace
%eigenfunctions of M.
%
%M = COMPUTE_EXTBASIS(M, k, order, normalize) Compute the extended basis of
%mesh M, conntaining the polynomials of the given order of the first k
%Laplace eigenfunctions of M and, if required, normalize it.

    % Check if M has an already computed basis
    if ~isfield(M, 'basis')
        % If k is not given, error
        if nargin < 2
            error(strcat("If M has not a precomputed Laplace basis, ", ...
                         "then a number of eigenfunctions if required."));
        end
        M = mesh.proc.compute_laplace_basis(M, k);
    end
    
    % If k is not given, then use all the eigenfunctions
    if nargin < 2
        k = size(M.basis, 2);
    end
    % If order is not given, then use 2
    if nargin < 3
        order = 2;
    end
    % If normalization is not given, don't do it
    if nargin < 4
        normalize = false;
    end
    
    % Compute the extended basis
    basis = M.basis(:, 1:k);
    M.extbasis = M.evecs(:, 1:k+1);
    prods = basis;
    for o = 2:order
        prods = repmat(basis, 1, k^(o - 1)) .* ...
                reshape(repmat(prods, k, 1), M.n, k^o);
        M.extbasis = [M.extbasis prods];
    end
    
    % If normalization is required, normalize the extended basis functions
    if ~normalize
        return;
    end
    % Memory-efficient normalization
    if size(M.extbasis, 2) < M.n
        norm_coeff = diag(M.extbasis' * M.A * M.extbasis);
    else
        norm_coeff = sum(M.extbasis' .* (M.A * M.extbasis)', 2);
    end
    M.extbasis = M.extbasis ./ repmat(sqrt(norm_coeff)', M.n, 1);
end

