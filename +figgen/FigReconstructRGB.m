% Plots the visual recostruction and the numerical error for a RGB signal
% on a shape using the following approaches:
%   - K Laplacian eigenfunctions
%   - N * K Laplacian eigenfunctions
%   - Eigenproducts
%   - Orthogonalized eigenproducts with high threshold
%   - Orthogonalized eigenproducts with low threshold and reiteration.
%
% This script reproduces the first half of Fig. 16 of the reference paper.
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
NumEigs = 200;
Ks = 20:20:100;
Order = 2;
OrthoThreshSmall = 1e-9;
OrthoThreshLarge = 1e-4;
FontSize = 48;


%% Load the mesh
fprintf("Loading mesh and computing eigenfunctions... ");
tic;
M = mesh.init(MeshFile);
M = mesh.transform.rotate(M, 'x', 90);
M = mesh.transform.rotate(M, 'z', -10);
[M.S, ~, M.A] = mesh.proc.FEM_higher(M, 1, 'Dirichlet');
[M.Phi, M.Lambda] = eigs(M.S, M.A, NumEigs, -1e-5);
M.Lambda = diag(M.Lambda);
F = load('texture_bunny.mat').M_col;
toc;


%% Plot
fig = figure;
fig.WindowState = 'maximized';

% Embed and reconstruct with eigs, prods and ortho
FK = cell(length(Ks), 1);
FP = cell(length(Ks), 1);
FS = cell(length(Ks), 1);
FL = cell(length(Ks), 1);

Ks = Ks(end);
for i = 1:length(Ks)
    fprintf("Iteration %d...", i);
    tic;
    
    % Plot original
    subplot(length(Ks), 5, (i - 1) * 5 + 1);
    utils.plot_rgb_map(M, F, true);
    
    % Eigs
    K = Ks(i);
    FK{i} = M.Phi(:, 1:K) * (M.Phi(:, 1:K)' * (M.A * F));
    subplot(length(Ks), 5, (i - 1) * 5 + 2);
    utils.plot_rgb_map(M, FK{i}, true);
    
    % Prods
    M.PolyPhi = eigprods(M, K - 1, Order, false);
    FP{i} = M.PolyPhi * (pinv(sqrt(M.A) * M.PolyPhi) * (sqrt(M.A) * F));
    subplot(length(Ks), 5, (i - 1) * 5 + 3);
    utils.plot_rgb_map(M, FP{i}, true);
    
    [M.Q, M.R] = orthogonalize(M.A, M.PolyPhi, OrthoThreshSmall, true);
    FS{i} = M.Q * (M.Q' * (M.A * F));
    subplot(length(Ks), 5, (i - 1) * 5 + 5);
    utils.plot_rgb_map(M, FS{i}, true);
    
    [M.Q, M.R] = orthogonalize(M.A, M.PolyPhi, OrthoThreshLarge, false);
    FL{i} = M.Q * (M.Q' * (M.A * F));
    subplot(length(Ks), 5, (i - 1) * 5 + 4);
    utils.plot_rgb_map(M, FL{i}, true);
    
    toc;
end














































































