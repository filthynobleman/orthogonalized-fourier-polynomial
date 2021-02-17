function fig = plot_scalar_function(M, func, cm, fig)
%PLOT_SCALAR_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

    T = M.TRIV;
    X = M.VERT(:, 1);
    Y = M.VERT(:, 2);
    Z = M.VERT(:, 3);
    if nargin < 4
        fig = figure;
    end
    if nargin < 3
        cm = jet;
    end
    if nargin < 2
        func = Z;
    end
    
    
    trisurf(T, X, Y, Z, func);
    colormap(cm);
end

