function [N, idxs] = cut_mesh(M, p, thresh, keepall, tokeep)
%CUT_MESH Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 4
        keepall = false;
    end

    if nargin < 5
        M.dijk = mesh.metrics.dijkstra(M, p);
        tokeep = find(M.dijk <= thresh);
    end
    
    if islogical(tokeep)
        tokeep = find(tokeep);
    end
    
    N.TRIV = M.TRIV((all(ismember(M.TRIV, tokeep), 2)), :);
    N.m = size(N.TRIV, 1);
    
    N.VERT = M.VERT(tokeep, :);
    N.X = N.VERT(:, 1);
    N.Y = N.VERT(:, 2);
    N.Z = N.VERT(:, 3);
    if isfield(M, 'dijk')
        N.dijk = M.dijk(tokeep);
    end
    N.n = size(N.VERT, 1);
    
    p = zeros(N.n, 1);
    p(tokeep) = 1:N.n;
    
    N.TRIV = p(N.TRIV);
    
    if keepall
        if isfield(M, 'A')
            N.A = M.A(:, tokeep);
            N.A = N.A(tokeep, :);
        end
        if isfield(M, 'S')
            N.S = M.S(:, tokeep);
            N.S = N.S(tokeep, :);
        end
        if isfield(M, 'basis')
            N.basis = M.basis(tokeep, :);
            N.basis = N.basis ./ repmat(sqrt(diag(N.basis' * N.A * N.basis))', N.n, 1);
        end
    end
    
    if nargout > 1
        idxs = tokeep;
    end
end

