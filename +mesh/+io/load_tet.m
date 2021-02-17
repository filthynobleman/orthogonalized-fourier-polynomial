function M = load_tet(filename)

    node = strcat(filename, ".node");
    ele = strcat(filename, ".ele");
    
    % Read the vertices
    fid = fopen(node);
    % First line contains informations
    vals = fscanf(fid, "%d  %d  %d  %d\n", [4 1]);
    M.Vn = vals(1);
    dim = vals(2);
    attr = vals(3);
    prop = vals(4);
    % Build the format properly
    format = "%d    ";
    for i = 1:(dim + attr + prop)
        format = strcat(format, "%f");
        if i < (dim + attr + prop)
            format = strcat(format, "  ");
        else
            format = strcat(format, "\n");
        end
    end
    % Read the vertices and exclude redundant informations
    M.V = fscanf(fid, format, [(1 + dim + attr + prop) M.Vn]);
    M.V = M.V(2:4, :)';
    fclose(fid);
    
    % Read the tetrahedrons
    fid = fopen(ele);
    % First line contains header
    vals = fscanf(fid, "%d  %d  %d\n", [3 1]);
    M.Tn = vals(1);
    dim = vals(2);
    attr = vals(3);
    % Build the format properly
    format = "%d    ";
    for i = 1:(dim + attr)
        format = strcat(format, "%d");
        if i < (dim + attr + prop)
            format = strcat(format, "  ");
        else
            format = strcat(format, "\n");
        end
    end
    % Read the tetrahedrons and exclude redundant informations
    M.T = fscanf(fid, format, [(1 + dim + attr) M.Tn]);
    M.T = uint32(M.T(2:5, :)') + 1;
    fclose(fid);
end

