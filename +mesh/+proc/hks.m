function HKS = hks(M, T, k)
%HKS Summary of this function goes here
%   Detailed explanation goes here


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
    
    
    HKS = zeros(M.n, length(T));
    for i = 1:length(T)
        HKS(:, i) = diag(mesh.proc.heat_kernel(M, T(i), k));
    end


end

