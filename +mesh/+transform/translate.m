function M = translate(M, translation)
%M = TRANSLATE(M, translation)
%   Translate the mesh in the 3-D space by the given 3-D vector.

    M.VERT = M.VERT + translation;
end

