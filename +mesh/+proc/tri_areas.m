function areas = tri_areas(M)
%areas = TRI_AREAS(M) This function computes the area of each triangle in
%the mesh and returns them in a column array. The value of the i-th cell is
%the area of the i-th triangle of M.
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 14 February 2020
    
    % Compute the edges' lengths
    E12 = sqrt(sum((M.VERT(M.TRIV(:, 1), :) - M.VERT(M.TRIV(:, 2), :)).^2, 2));
    E23 = sqrt(sum((M.VERT(M.TRIV(:, 2), :) - M.VERT(M.TRIV(:, 3), :)).^2, 2));
    E13 = sqrt(sum((M.VERT(M.TRIV(:, 1), :) - M.VERT(M.TRIV(:, 3), :)).^2, 2));
    % Compute the semiperemetrs
    S = (E12 + E23 + E13) / 2;
    % Compute the areas
    areas = sqrt(S .* (S - E12) .* (S - E23) .* (S - E13));
end

