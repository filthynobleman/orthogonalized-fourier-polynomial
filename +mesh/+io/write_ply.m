function write_ply(filename, M, func, cm, center_cmap, format)
%WRITE_PLY(filename, M) Writes the mesh M to the given filename in binary
%PLY file format.
%
%WRITE_PLY(filename, M, func, cm) Also associates the mesh with a vertex
%coloring determined by the scalar function and the colormap.
%
%WRITE_PLY(filename, M, func, cm, center_cmap) Also determines if the
%colormap must be centered at the zero value of the function. Default is
%true.
%
%WRITE_PLY(filename, M, func, cm, center_cmap, format) Also allows to
%choose the format, between ASCII and binary.
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 17 February 2020

    if nargin < 5
        center_cmap = true;
    end
    if nargin < 6
        format = "binary";
    end

    if lower(format) == "ascii"
        write_ply_ascii(filename, M, func, cm, center_cmap);
    else
        if nargin < 3
            write_ply_binary(filename, M);
        else
            write_ply_binary(filename, M, func, cm, center_cmap);
        end
    end
end



function write_ply_ascii(filename, M, func, cm, center_cmap)
%WRITE_PLY(filename, M) Writes the mesh M to the given filename in ASCII
%PLY file format.
%
%WRITE_PLY(filename, M, func, cm) Also associates the mesh with a vertex
%coloring determined by the scalar function and the colormap.
%
%WRITE_PLY(filename, M, func, cm, center_cmap) Also determines if the
%colormap must be centered at the zero value of the function. Default is
%true.
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 17 February 2020

    fid = fopen(filename, "w");
    if fid == 0
        error("Cannot open file %s.", filename);
    end

    % Write the header
    fwrite(fid, "ply", 'char');
    fwrite(fid, newline);
    fwrite(fid, "format ascii 1.0", 'char');
    fwrite(fid, newline);
    % Coordinates
    fwrite(fid, sprintf("element vertex %d", M.n));
    fwrite(fid, newline);
    fwrite(fid, "property float x", 'char');
    fwrite(fid, newline);
    fwrite(fid, "property float y", 'char');
    fwrite(fid, newline);
    fwrite(fid, "property float z", 'char');
    fwrite(fid, newline);
    % If a vertex coloring has been given, add the properties
    if nargin > 2
        fwrite(fid, "property uchar red", 'char');
        fwrite(fid, newline);
        fwrite(fid, "property uchar green", 'char');
        fwrite(fid, newline);
        fwrite(fid, "property uchar blue", 'char');
        fwrite(fid, newline);
    end
    % Triangles
    fwrite(fid, sprintf("element face %d", M.m), 'char');
    fwrite(fid, newline);
    fwrite(fid, "property list uchar int vertex_indices", 'char');
    fwrite(fid, newline);
    fwrite(fid, "end_header", 'char');
    fwrite(fid, newline);
    
    if nargin > 2
        % If a coloring has been provided, add vertices and colors
        f = func;
        if center_cmap
            f = f + max(abs(func));
        end
        f = f / (max(f) - min(f));
        f = f * (size(cm, 1) - 1) + 1;
        RGB = cm(round(f), :);
        RGB = RGB * 255;
        fprintf(fid, "%.6f %.6f %.6f %.0f %.0f %.0f\n", [M.VERT RGB]');
    else
        % Otherwise add only the vertices
        fprintf(fid, "%.6f %.6f %.6f\n", M.VERT');
    end
    
    % Write the triangles
    fprintf(fid, "3 %.0f %.0f %.0f\n", M.TRIV' - 1);
    fclose(fid);
end


function write_ply_binary(filename, M, func, cm, center_cmap)
%WRITE_PLY(filename, M) Writes the mesh M to the given filename in binary
%PLY file format.
%
%WRITE_PLY(filename, M, func, cm) Also associates the mesh with a vertex
%coloring determined by the scalar function and the colormap.
%
%WRITE_PLY(filename, M, func, cm, center_cmap) Also determines if the
%colormap must be centered at the zero value of the function. Default is
%true.
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 18 February 2020

    fid = fopen(filename, "wb");
    if fid == 0
        error("Cannot open file %s.", filename);
    end

    % Write the header
    fwrite(fid, "ply", 'char');
    fwrite(fid, newline, 'char');
    fwrite(fid, "format binary_little_endian 1.0", 'char');
    fwrite(fid, newline, 'char');
    % Coordinates
    fwrite(fid, sprintf("element vertex %d", M.n), 'char');
    fwrite(fid, newline, 'char');
    fwrite(fid, "property float x", 'char');
    fwrite(fid, newline, 'char');
    fwrite(fid, "property float y", 'char');
    fwrite(fid, newline, 'char');
    fwrite(fid, "property float z", 'char');
    fwrite(fid, newline, 'char');
    % If a vertex coloring has been given, add the properties
    if nargin > 2
        fwrite(fid, "property uchar red", 'char');
        fwrite(fid, newline, 'char');
        fwrite(fid, "property uchar green", 'char');
        fwrite(fid, newline, 'char');
        fwrite(fid, "property uchar blue", 'char');
        fwrite(fid, newline, 'char');
    end
    % Triangles
    fwrite(fid, sprintf("element face %d", M.m), 'char');
    fwrite(fid, newline, 'char');
    fwrite(fid, "property list uchar int vertex_indices", 'char');
    fwrite(fid, newline, 'char');
    fwrite(fid, "end_header", 'char');
    fwrite(fid, newline, 'char');
    
%     fclose(fid);
    
    M.VERT = single(M.VERT');
    M.TRIV = int32(M.TRIV' - 1);
    if nargin > 2
        % If a coloring has been provided, add vertices and colors
        f = func;
        if center_cmap
            f = [f; max(abs(f)); -max(abs(f))];
%             f = f + max(abs(func));
        end
        f = f - min(f);
        f = f / max(f);
        f = f * (size(cm, 1) - 1) + 1;
        if center_cmap
            f = f(1:end - 2);
        end
        RGB = cm(round(f), :);
        RGB = uint8(RGB * 255)';
%         mesh.io.write_ply_verts_rgb(M.VERT', RGB', filename);
        buf = mesh.io.create_ply_content(M.VERT, M.TRIV, RGB);
    else
        % Otherwise add only the vertices
%         fwrite(fid, M.VERT', 'float32');
        buf = mesh.io.create_ply_content(M.VERT, M.TRIV);
    end
    
    fwrite(fid, buf, 'uint8');
    fclose(fid);
    
    % Write the triangles
%     mesh.io.write_ply_triv(int32(M.TRIV'), filename);
%     fclose(fid);
end