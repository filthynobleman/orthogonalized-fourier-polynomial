function M = scale(M, factor)
%M = SCALE(M, factor)
%   Scales the vertices coordinates by the given factor.

    M.VERT = M.VERT * factor;
    M.X = M.VERT(:, 1);
    M.Y = M.VERT(:, 2);
    M.Z = M.VERT(:, 3);
end

