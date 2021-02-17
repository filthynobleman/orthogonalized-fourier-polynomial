function [N, idxs] = get_area(M, Start, Area)

    E = [M.TRIV(:, [1, 2]); M.TRIV(:, [2, 3]); M.TRIV(:, [3, 1])];
    Adj = sparse(E(:, 1), E(:, 2), 1, M.n, M.n);
    Adj = double(Adj + Adj' > 0);
    Adj(1:M.n+1:end) = 1;

    f = zeros(M.n, 1);
    f(Start) = 1;
    
    N.Area = 0;
    while N.Area < Area
        f = Adj * f;
        N = mesh.proc.cut_mesh(M, 1, 1, false, f > 0);
        N.Area = sum(mesh.proc.tri_areas(N));
    end
    
    if nargout > 1
        idxs = find(f > 0);
    end
end

