function M = load_off(filename)
%M = LOAD_OFF(filename) Loads a mesh from the given OFF file.
%
%   The file is expected to begin directly after the header, without
%   containing black or commented lines.
%   The function loads a triangle mesh, so the file is expected to contain
%   only triangles.

    % Ensure the function got a string in input
    if ~isstring(filename) && ~ischar(filename)
        errmsg = "Function load_off() only gets ";
        errmsg = strcat(errmsg, "strings and character arrays.");
        error(errmsg);
    end
    
    % Initialize the mesh
    M = [];
    
    % Try to open the file and check if it has benn successfully opened
    fid = fopen(filename, 'r');
    if fid == -1
        error("Cannot open " + filename + " for reading.");
    end
    
    % Check the file header to ensure the given file is an OFF
    if fscanf(fid, "%s\n", 1) ~= "OFF"
        error("File " + filename + " is not a valid OFF file.");
    end
    
    % Ignore #-starting and empty lines
    while true
        line = fgets(fid);
        if startsWith(line, "#")
            continue;
        elseif length(line) < 3
            continue;
        else
            break;
        end
    end
    
    % Read the header
    header = sscanf(line, "%d %d %d\n", [3 1]);
    M.n = header(1);
    M.m = header(2);
    % Read the vertices and the triangles
    M.VERT = fscanf(fid, "%f %f %f\n", [3 M.n])';
    M.TRIV = fscanf(fid, "3 %d %d %d\n", [3 M.m])';
    if any(M.TRIV == 0, 'all')
        M.TRIV = M.TRIV + 1;
    end
    
    % Close the file
    fclose(fid);
end

