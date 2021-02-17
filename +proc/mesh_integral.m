function varargout = mesh_integral(M, f)
%fint = MESH_INTEGRAL(M, f) Computes the integral of the function f over
%the mesh. If f is a matrix, it is treated as a set of functions, one for
%each column.
%
%[fint, ftris] = MESH_INTEGRAL(M, f) the same as before, but it also
%returns the value of f on each triangle.

    % Compute the areas of the triangles
    areas = mesh.proc.tri_areas(M);
    % Compute the functions at the vertices
    Fx1 = f(M.TRIV(:, 1), :);
    Fx2 = f(M.TRIV(:, 2), :);
    Fx3 = f(M.TRIV(:, 3), :);
    % Compute the functions over the single triangles
    ftris = (Fx1 + Fx2 + Fx3) .* areas / 3;
    % Compute the integrals over the mesh
    fint = sum(ftris);
    
    varargout{1} = fint;
    if nargout == 2
        varargout{2} = ftris;
    end

end

