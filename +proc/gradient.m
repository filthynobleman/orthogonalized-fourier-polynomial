function grad = gradient(M, f, onverts)
%grad = GRADIENT(M, f) 
%Computes the gradient of function f on the triangles of mesh M.
%
%grad = GRADIENT(M, f, onverts) If onverts is true, then computes the
%gradient of f on vertices, rather than on triangles.

    % Compute the lengths of the edges
    E12 = M.VERT(M.TRIV(:, 2), :) - M.VERT(M.TRIV(:, 1), :);
    E13 = M.VERT(M.TRIV(:, 3), :) - M.VERT(M.TRIV(:, 1), :);
    % Compute the gradient of f tilde
    F12 = f(M.TRIV(:, 2)) - f(M.TRIV(:, 1));
    F13 = f(M.TRIV(:, 3)) - f(M.TRIV(:, 1));

    % Compute the gradient with the MEX function
    grad = proc.compute_gradient(E12, E13, F12, F13, uint32(M.m));
    
    if nargin > 2 && onverts
        grad = proc.grad_on_verts(M, grad);
    end
end

