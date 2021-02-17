% Plots the first eigenfunctions, eigenproducts and orthogonalized
% eigenproducts on a given shape.
%
% This script reproduces Fig. 8 of the reference paper.
%
% This code is part of the implementation of the paper "Orthogonalized 
% Fourier Polynomials for Signal Approximation and Transfer", authored by
% Filippo Maggioli, Simone Melzi, Michael Bronstein, Maks Ovsjanikov and
% Emanuele Rodolà.

%% Clear the workspace
clear all
close all
clc


%% Parameters
MeshFile = 'Meshes/bunny.off';
NumRows = 1;
NumCols = 5;
OrthoThresh = 1e-9;
FontSize = 48;


%% Load the mesh
fprintf("Loading mesh and computing eigenfunctions... ");
tic;
M = mesh.init(MeshFile);
[M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Dirichlet');
NumFuncs = NumRows * NumCols;
[M.Phi, M.Lambda] = eigs(M.S, M.A, NumFuncs, -1e-5);
M.Lambda = diag(M.Lambda);
toc;


%% Compute eigenproducts
fprintf("Computing eigenproducts... ");
tic;
M.PolyPhi = eigprods(M, NumFuncs - 1, 2, false);
toc;


%% Compute orthogonal basis
fprintf("Orthogonalization... ");
tic;
[M.Q, M.R] = orthogonalize(M.A, M.PolyPhi, OrthoThresh, true);
toc;



%% Plot
% Plot eigenfunctions
for i = 1:NumRows
    for j = 1:NumCols
        Idx = (i - 1) * NumCols + j;
        PlotIdx = Idx;
        subplot(3 * NumRows, NumCols, PlotIdx);
        utils.plot_scalar_map(M, M.Phi(:, Idx), utils.cmaps.bwr, true);
        title(sprintf("$\\varphi_{%d}$", Idx - 1), ...
              'Interpreter', 'LaTeX', ...
              'FontSize', FontSize);
    end
end

% Plot eigenproducts
for i = 1:NumRows
    for j = 1:NumCols
        Idx = (i - 1) * NumCols + j;
        PlotIdx = NumRows * NumCols + Idx;
        subplot(3 * NumRows, NumCols, PlotIdx);
        utils.plot_scalar_map(M, M.PolyPhi(:, NumFuncs + Idx), utils.cmaps.bwr, true);
        Idx1 = floor(Idx / NumFuncs) + 1;
        Idx2 = mod(Idx, NumFuncs) + 2 * (Idx1 - 1);
        title(sprintf("$\\varphi_{%d}\\varphi_{%d}$", Idx1, Idx2), ...
              'Interpreter', 'LaTeX', ...
              'FontSize', FontSize);
    end
end

% Plot orthogonal
for i = 1:NumRows
    for j = 1:NumCols
        Idx = (i - 1) * NumCols + j;
        PlotIdx = 2 * NumRows * NumCols + Idx;
        subplot(3 * NumRows, NumCols, PlotIdx);
        utils.plot_scalar_map(M, M.Q(:, NumFuncs + Idx), utils.cmaps.bwr, true);
        title(sprintf("$Q_{%d}$", Idx - 1), ...
              'Interpreter', 'LaTeX', ...
              'FontSize', FontSize);
    end
end














































































