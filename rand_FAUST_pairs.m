function MeshFiles = rand_FAUST_pairs(indir)
%RAND_FAUST Summary of this function goes here
%   Detailed explanation goes here


    if nargin == 0
        indir = "Meshes";
    end
    
    L = 20;
    MeshFiles = cell(L, 2);
    for i = 1:L
        One.Shape = randi([0, 9]);
        One.Pose = randi([0, 9]);
        One.MeshFile = sprintf("%s/tr_reg_0%d%d.ply", indir, One.Shape, One.Pose);
        Two.Shape = randi([0, 9]);
        Two.Pose = randi([0, 9]);
        Two.MeshFile = sprintf("%s/tr_reg_0%d%d.ply", indir, Two.Shape, Two.Pose);
        MeshFiles{i, 1} = One.MeshFile;
        MeshFiles{i, 2} = Two.MeshFile;
    end
end

