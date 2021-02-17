% Reconstruction of the XYZ coordinates function for two shapes, using a
% linear combination of eigenproducts with increasing order and a limited
% number of eigenfunctions.
%
% This script reproduces Fig. 1 of the reference paper.
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
NumEigs = 5;
OrthoOrder = 1:5;
SVDThresh = (NumEigs - 2) * 2e-3;
MeshFiles = ["Meshes/189_filigree.off";
             "Meshes/bimba.off"];
FontSize = 24;
ConstCol = 0.95;
Colors = [0 0 0;
          ConstCol 0 0;
          ConstCol 0.5 0;
          ConstCol ConstCol 0.5;
          ConstCol ConstCol ConstCol];
ColMap = utils.cmaps.blend(flipud(Colors), 256);
ConstCM = utils.cmaps.constant(ConstCol);


%% Load the mesh
fprintf("laoding the mesh... ");
tic;
M = mesh.init(MeshFiles(1));
M = mesh.transform.rotate(M, 'x', 90);
M = mesh.transform.rotate(M, 'z', -60);
toc;


%% Solve the eigenproblem
fprintf("Computing stiff/mass matrices... ");
tic;
[M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Dirichlet');
toc;

fprintf("Solving the eigenproblem for K = %d... ", NumEigs);
tic;
[M.Phi, M.Lambda] = eigs(M.S, M.A, NumEigs + 1, -1e-5);
M.Lambda = diag(M.Lambda);
toc;


%% Compute eigenproducts and orthogonal basis

fprintf("Computing orthogonal bases...\n");
M.Q = cell(length(OrthoOrder), 1);
for i = 1:length(OrthoOrder)
    N = OrthoOrder(i);
    fprintf("Order %d... ", N);
    tic;
    Tmp = eigprods(M, NumEigs, N, false);
    M.Q{i} = orthogonalize(M.A, Tmp, 1e-9, true);
    toc;
end
clear Tmp;



%% Plot the original shape
fig = figure;
fig.WindowState = 'maximized';


subplot(2, 6, 1);
utils.plot_scalar_map(M, zeros(M.n, 1), ConstCM, true, 'interp', false);
title("Original", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize);
  
  
%% Reconstruct with orthogonal
Lim = 0;
for i = 1:length(OrthoOrder)
    Emb = M.Q{i}' * (M.A * M.VERT);
    VertRec = M.Q{i} * Emb;
    
    N = M;
    N.VERT = VertRec;
    N.X = N.VERT(:, 1);
    N.Y = N.VERT(:, 2);
    N.Z = N.VERT(:, 3);
    subplot(2, 6, i + 1);
    Den = vecnorm(M.VERT, 2, 2);
    Err = vecnorm(N.VERT - M.VERT, 2, 2) ./ Den;
    Lim = max(Lim, max(abs(Err)));
    utils.plot_scalar_map(N, Err, ColMap, true, 'interp', false);
    TitleString = {sprintf("Orthogonal Order %d", OrthoOrder(i));
                   sprintf("Err: %.3f", sqrt((Err' * M.A * Err) / (Den' * M.A * Den)))};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


for i = 1:length(OrthoOrder)
    subplot(2, 6, i + 1);
    caxis([0, Lim]);
end







%% Load the mesh
fprintf("laoding the mesh... ");
tic;
M = mesh.init(MeshFiles(2));
M = mesh.transform.rotate(M, 'x', 80);
M = mesh.transform.rotate(M, 'z', -45);
toc;


%% Solve the eigenproblem
fprintf("Computing stiff/mass matrices... ");
tic;
[M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Dirichlet');
toc;

fprintf("Solving the eigenproblem for K = %d... ", NumEigs);
tic;
[M.Phi, M.Lambda] = eigs(M.S, M.A, NumEigs + 1, -1e-5);
M.Lambda = diag(M.Lambda);
toc;


%% Compute eigenproducts and orthogonal basis

fprintf("Computing orthogonal bases...\n");
M.Q = cell(length(OrthoOrder), 1);
for i = 1:length(OrthoOrder)
    N = OrthoOrder(i);
    fprintf("Order %d... ", N);
    tic;
    Tmp = eigprods(M, NumEigs, N, false);
    M.Q{i} = orthogonalize(M.A, Tmp, 1e-9, true);
    toc;
end
clear Tmp;



%% Plot the original shape
subplot(2, 6, 7);
utils.plot_scalar_map(M, zeros(M.n, 1), ConstCM, true, 'interp', false);
title("Original", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize);
  
  
%% Reconstruct with orthogonal
Lim = 0;
for i = 1:length(OrthoOrder)
    Emb = M.Q{i}' * (M.A * M.VERT);
    VertRec = M.Q{i} * Emb;
    
    N = M;
    N.VERT = VertRec;
    N.X = N.VERT(:, 1);
    N.Y = N.VERT(:, 2);
    N.Z = N.VERT(:, 3);
    subplot(2, 6, 6 + i + 1);
    Den = vecnorm(M.VERT, 2, 2);
    Err = vecnorm(N.VERT - M.VERT, 2, 2) ./ Den;
    Lim = max(Lim, max(abs(Err)));
    utils.plot_scalar_map(N, Err, ColMap, true, 'interp', false);
    TitleString = {sprintf("Orthogonal Order %d", OrthoOrder(i));
                   sprintf("Err: %.3f", sqrt((Err' * M.A * Err) / (Den' * M.A * Den)))};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


for i = 1:length(OrthoOrder)
    subplot(2, 6, 6 + i + 1);
    caxis([0, Lim]);
end






















































