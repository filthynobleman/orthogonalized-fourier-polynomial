function PolyPhi = eigprods(M, k, n, normalized)
%PolyPhi = EIGPRODS(M, k) Computes the second order polynomial basis from
%the first non-constant k Laplacian eigenfunctions of mesh M.
%
%PolyPhi = EIGPRODS(M, k, n) Computes the n-th ordcer polynomial basis.
%
%PolyPhi = EIGPRODS(M, k, n, normalized) Decide if normalize or not the
%polynomial basis so that the functions have squared norm equals to 1 on
%the mesh. Default is false.
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 11 September 2020

    if nargin < 4
        normalized = false;
    end
    if nargin < 3
        n = 2;
    end
    
    % Check for the mass and stiffness matrices
    if ~isfield(M, 'A') && ~isfield(M, 'S')
        [M.S, M.A] = mesh.proc.FEM_higher(M, 1, 'Dirichlet');
    end
    
    % Check if the basis has already been computed and has enough elements
    if ~isfield(M, 'Phi') || size(M.Phi, 2) < k + 1
        [M.Phi, M.Lambda] = eigs(M.S, M.A, k + 1, -1e-5);
        M.Lambda = diag(M.Lambda);
    end
    
    PolyPhi = M.Phi(:, 1:k+1);
    Basis = PolyPhi(:, 2:end);
    Prods = Basis;
    for i = 2:n
        tmp1 = repmat(Basis, 1, k^(i - 1));
        tmp2 = reshape(repmat(Prods, k, 1), size(M.A, 1), k^i);
        Prods = (tmp1 .* tmp2);
        PolyPhi = [PolyPhi, Prods];
    end
    PolyPhi = PolyPhi(:, sub_index(n, k));
    
    
    if normalized
        norms = diag(PolyPhi' * M.A * PolyPhi);
        PolyPhi = PolyPhi ./ sqrt(norms');
    end
end





function [subs_order] = sub_index(order, k)

    PolyPhi = 0:k;
    Basis = 1:k;
    Prods = Basis;
    for i = 2:order
        tmp1 = repmat(Basis, 1, k^(i - 1));
        tmp2 = reshape(repmat(Prods, k, 1), i - 1, k^i);
        Prods = [tmp1; tmp2];
        PolyPhi = [PolyPhi; zeros(1, length(PolyPhi))];
        PolyPhi = [PolyPhi, Prods];
    end
    PolyPhi = sort(PolyPhi, 1);
    [~, subs_order, ~] = unique(PolyPhi', 'rows');
end