function write_off(M, filename, fp_prec)
%WRITE_OFF(M, filename, fp_prec) This function saves the given mesh as an
%OFF file using the given floating point precision for representing double
%values.
%
%   This function gets in input a mesh a string, representing a valid
%   filename to a writable position, and an integer. It saves the given
%   mesh as an OFF file to the given filename.
%   The double values representing the vertices' coordinates are saved
%   using a floating point precision defined by the given integer.
%
%
%WRITE_OFF(M, filename) This function saves the given mesh as an OFF file.
%
%   This function gets in input a mesh and a string, representing a valid
%   filename to a writable position. It saves the given mesh as an OFF file
%   to the given filename.
%   It is equivalent to call write_off(M, path, 6).
    
    % Ensure the function got a string in input
    if ~isstring(filename) && ~ischar(filename)
        errmsg = "Function write_off() only gets ";
        errmsg = strcat(errmsg, "strings and ...character arrays.");
        error(errmsg);
    end
    
    % Check for the third argument. If it is not given, the default
    % floating point precision for the representation is 6
    if nargin < 3
        fp_prec = 6;
    end
    
    % Try to open the file
    fid = fopen(filename, 'w');
    if fid == -1
        error("Cannot open " + filename + " for writing.");
    end
    % Write the header
    fprintf(fid, "OFF\n");
    fprintf(fid, "%d %d 0\n", M.n, M.m);
    % Write the vertices
    formatSpec = "%." + int2str(fp_prec) + "f";
    formatSpec = formatSpec + ' ' + formatSpec + ' ' + formatSpec + '\n';
    fprintf(fid, formatSpec, M.VERT);
    % Write the triangles
    fprintf(fid, "3 %d %d %d\n", M.TRIV);
    % Close the file
    fclose(fid);
end

