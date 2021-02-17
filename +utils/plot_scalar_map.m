function fig = plot_scalar_map(M, f, cm, shadow, shading_opt, center_map)
%PLOT_SCALAR_MAP(M) Plots the mesh M in white.
%
%PLOT_SCALAR_MAP(M, f, cm) Plots the scalar map f over the mesh M, using
%the colormap cm.
%Default values for f is an all zero vector.
%Default value for cm is white.
%
%PLOT_SCALAR_MAP(M, f, cm, shadow) Also plots a shadow of the mesh on the
%floor.
%Default value is false.
%
%PLOT_SCALAR_MAP(M, f, cm, shadow, shading_opt) Also allows to define the
%shading options for the plot. For valid values check 'help shading'.
%Default value is 'interp'.
%
%PLOT_SCALAR_MAP(M, f, cm, shadow, shading_opt, center_map) Also centers
%the colormap so that the middle value of the map points at zero.
%Default value is true.
%
%fig = PLOT_SCALAR_MAP(--) The function creates a new figure handle, plots
%the mesh and returns the handle.
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 14 February 2020

    if nargout > 0
        fig = figure;
    end
    
    if nargin == 1
        f = zeros(M.n, 1);
        cm = white;
    end
    if nargin < 4
        shadow = false;
    end
    if nargin < 5
        shading_opt = 'interp';
    end
    if nargin < 6
        center_map = true;
    end
    
    if ~isfield(M, 'X')
        M.X = M.VERT(:, 1);
    end
    if ~isfield(M, 'Y')
        M.Y = M.VERT(:, 2);
    end
    if ~isfield(M, 'Z')
        M.Z = M.VERT(:, 3);
    end
    
    if ~shadow
        trisurf(M.TRIV, M.X, M.Y, M.Z, f, ...
            'SpecularStrength', 0.1, ...
            'DiffuseStrength', 0.35);
        axis off;
        axis equal;
        shading(shading_opt);
        colormap(gca, cm);
        view([0, 90]);
        light;
        lighting phong;
        camlight head;
    else
%         M = mesh.transform.rotate(M, 'x', 90);
        M = mesh.transform.rotate(M, 'z', 30);
        addpath('func_render/');
        M.name = 'mesh';
        M.nv = M.n;
        M.surface = M;
%         MinF = min(f);
%         MaxF = max(f);
%         NormF = (f - MinF) / (MaxF - MinF);
%         FIdx = min(round(NormF .* (length(cm) - 1)) + 1, length(cm));
        render_mesh(M,...
            'MeshVtxColor', f,... % [R G B] for each of the vertex
            'CameraPos', [0,90],...
            'FaceAlpha', 0.9);
        shading(shading_opt);
        colormap(gca, cm);
        view([0, 10]);
    end
    if center_map
        lim = max(abs(f));
        caxis([-lim, lim]);
    end
end

