function gradv = grad_on_verts(M, grad)
%gradv = GRAD_ON_VERTS(M, grad)
%Given a gradient grad on the triangles of the mesh M, compute the same
%gradient, applied to the vertices of M.

    areas = mesh.proc.tri_areas(M);
    P = sparse(M.TRIV(:), [1:M.m 1:M.m 1:M.m], repmat(areas, 3, 1));
    gradv = sparse(1:M.n, 1:M.n, 1 ./ sum(P, 2)) * (P * grad);
end

