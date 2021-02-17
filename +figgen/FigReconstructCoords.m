% Plots the visual recostruction and the numerical error for the XYZ
% coordinates functions on a shape using the following approaches:
%   - K Laplacian eigenfunctions
%   - N * K Laplacian eigenfunctions
%   - Eigenproducts
%   - Orthogonalized eigenproducts with high threshold
%   - Orthogonalized eigenproducts with low threshold and reiteration.
%
% This script reproduces the first half of Fig. 17 of the reference paper.
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
MeshFile = 'Meshes/205-Greek_Sculpture.off';
K = 50;
N = 2;
OrthoThresh = 1e-2;
OrthoThreshReiter = 1e-9;
SVDThresh = (K - 2) * 2e-3;
FontSize = 24;
PlotTitles = true;
XRotation = 80; %72.5;
YRotation = 0.0;
ZRotation = 0.0;


%% Load the meshe and compute the needed fields
fprintf("Loading the meshes... ");
tic;
M = mesh.init(MeshFile);
MPlot = mesh.transform.rotate(M, 'x', XRotation);
MPlot = mesh.transform.rotate(MPlot, 'y', YRotation);
MPlot = mesh.transform.rotate(MPlot, 'z', ZRotation);
% Src = mesh.transform.normalize(Src);
toc;

fprintf("Solve the eigenproblem... ");
tic;
[M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Dirichlet');
[M.Phi, M.Lambda] = eigs(M.S, M.A, N * K, -1e-5);
M.Lambda = diag(M.Lambda);
toc;


%% Initialize a function to transfer and plot it
f = zeros(M.n, 1);
Lim = max(abs(f));

subplot(2, 3, 1);
utils.plot_scalar_map(MPlot, f, white, true);
if PlotTitles
    title('Source', ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end

%% Compute NK functional map and transfer
Emb = M.Phi(:, 1:N*K)' * (M.A * M.VERT);
FRec = M.Phi(:, 1:N*K) * Emb;
Diff = vecnorm(M.VERT - FRec, 2, 2);
VertNorm = vecnorm(M.VERT, 2, 2);
ErrNK = sqrt(Diff' * M.A * Diff) / sqrt(VertNorm' * M.A * VertNorm);
Lim = max(Lim, max(abs(Diff)));

subplot(2, 3, 2);
New = M;
New.VERT = FRec;
New = mesh.transform.rotate(New, 'x', XRotation);
New = mesh.transform.rotate(New, 'y', YRotation);
New = mesh.transform.rotate(New, 'z', ZRotation);
utils.plot_scalar_map(New, Diff, flipud(hot), true, 'interp', false);
if PlotTitles
    TitleString = {"$NK$ Eigs";
                   sprintf("Err: %.3f", ErrNK)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Compute K functional map and transfer
Emb = M.Phi(:, 1:K)' * (M.A * M.VERT);
FRec = M.Phi(:, 1:K) * Emb;
Diff = vecnorm(M.VERT - FRec, 2, 2);
VertNorm = vecnorm(M.VERT, 2, 2);
ErrK = sqrt(Diff' * M.A * Diff) / sqrt(VertNorm' * M.A * VertNorm);
Lim = max(Lim, max(abs(Diff)));

subplot(2, 3, 3);
New = M;
New.VERT = FRec;
New = mesh.transform.rotate(New, 'x', XRotation);
New = mesh.transform.rotate(New, 'y', YRotation);
New = mesh.transform.rotate(New, 'z', ZRotation);
utils.plot_scalar_map(New, Diff, flipud(hot), true, 'interp', false);
if PlotTitles
    TitleString = {"$K$ Eigs";
                   sprintf("Err: %.3f", ErrK)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Transfer with eigenproducts
M.PolyPhi = eigprods(M, K - 1, N, false);

Emb = pinv(sqrt(M.A) * M.PolyPhi) * (sqrt(M.A) * M.VERT);
FRec = M.PolyPhi * Emb;
Diff = vecnorm(M.VERT - FRec, 2, 2);
VertNorm = vecnorm(M.VERT, 2, 2);
ErrPoly = sqrt(Diff' * M.A * Diff) / sqrt(VertNorm' * M.A * VertNorm);
Lim = max(Lim, max(abs(Diff)));

subplot(2, 3, 4);
New = M;
New.VERT = FRec;
New = mesh.transform.rotate(New, 'x', XRotation);
New = mesh.transform.rotate(New, 'y', YRotation);
New = mesh.transform.rotate(New, 'z', ZRotation);
utils.plot_scalar_map(New, Diff, flipud(hot), true, 'interp', false);
if PlotTitles
    TitleString = {"Poly";
                   sprintf("Err: %.3f", ErrPoly)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Transfer with orthogonal basis
M.Q = orthogonalize(M.A, M.PolyPhi, OrthoThresh, false);

Emb = M.Q' * (M.A * M.VERT);
FRec = M.Q * Emb;
Diff = vecnorm(M.VERT - FRec, 2, 2);
VertNorm = vecnorm(M.VERT, 2, 2);
ErrOrtho = sqrt(Diff' * M.A * Diff) / sqrt(VertNorm' * M.A * VertNorm);
Lim = max(Lim, max(abs(Diff)));

subplot(2, 3, 5);
New = M;
New.VERT = FRec;
New = mesh.transform.rotate(New, 'x', XRotation);
New = mesh.transform.rotate(New, 'y', YRotation);
New = mesh.transform.rotate(New, 'z', ZRotation);
utils.plot_scalar_map(New, Diff, flipud(hot), true, 'interp', false);
if PlotTitles
    TitleString = {"Ortho";
                   sprintf("Err: %.3f", ErrOrtho)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Transfer with stressed orthogonal basis
M.QS = orthogonalize(M.A, M.PolyPhi, OrthoThreshReiter, true);

Emb = M.QS' * (M.A * M.VERT);
FRec = M.QS * Emb;
Diff = vecnorm(M.VERT - FRec, 2, 2);
VertNorm = vecnorm(M.VERT, 2, 2);
ErrOrthoReiter = sqrt(Diff' * M.A * Diff) / sqrt(VertNorm' * M.A * VertNorm);
Lim = max(Lim, max(abs(Diff)));

subplot(2, 3, 6);
New = M;
New.VERT = FRec;
New = mesh.transform.rotate(New, 'x', XRotation);
New = mesh.transform.rotate(New, 'y', YRotation);
New = mesh.transform.rotate(New, 'z', ZRotation);
utils.plot_scalar_map(New, Diff, flipud(hot), true, 'interp', false);
if PlotTitles
    TitleString = {"Ortho Reiter";
                   sprintf("Err: %.3f", ErrOrthoReiter)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end



%% Caxis all
for i = 1:6
    subplot(2, 3, i);
    caxis([0, Lim]);
end






































































