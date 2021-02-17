function MeshFiles = rand_pairs(indir)
%RAND_FAUST Summary of this function goes here
%   Detailed explanation goes here


    if nargin == 0
        indir = "Meshes";
    end
    DirContent = dir(sprintf("%s/*.ply", indir));
    NumFiles = length(DirContent);
    
    L = 20;
    MeshFiles = cell(L, 2);
    for i = 1:L
        MeshFiles{i, 1} = sprintf("%s/%s", indir, DirContent(randi([1, NumFiles])).name);
        MeshFiles{i, 2} = sprintf("%s/%s", indir, DirContent(randi([1, NumFiles])).name);
    end
end

