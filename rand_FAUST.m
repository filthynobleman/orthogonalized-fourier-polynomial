function MeshFiles = rand_FAUST(indir)
%RAND_FAUST Summary of this function goes here
%   Detailed explanation goes here


    if nargin == 0
        indir = "Meshes";
    end
    
    MeshFiles = {};
    for i = 0:9
        idx = datasample(0:9, 2, 'Replace', false);
        MeshFiles{end + 1} = sprintf("%s/tr_reg_0%d%d.ply", indir, i, idx(1));
        MeshFiles{end + 1} = sprintf("%s/tr_reg_0%d%d.ply", indir, i, idx(2));
    end
    MeshFiles = MeshFiles';
end

