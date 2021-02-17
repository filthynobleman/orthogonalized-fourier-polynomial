function M = load_obj(filename)
%LOAD_OBJ Summary of this function goes here
%   Detailed explanation goes here


%     OBJ = mesh.io.read_wobj(filename);
%     M.VERT = OBJ.vertices;
%     M.TRIV = OBJ.objects(3).data.vertices;

    [M.VERT, M.TRIV, ~] = mesh.io.loadawobj(filename);
    if size(M.VERT, 2) > size(M.VERT, 1)
        M.VERT = M.VERT';
    end
    if size(M.TRIV, 2) > size(M.TRIV, 1)
        M.TRIV = M.TRIV';
    end
    M.n = size(M.VERT, 1);
    M.m = size(M.TRIV, 1);
end

