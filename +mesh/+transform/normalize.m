function M = normalize(M)
%NORMALIZE Summary of this function goes here
%   Detailed explanation goes here


    M = mesh.transform.scale(M, 1 / sqrt(sum(mesh.proc.tri_areas(M))));

end

