% Example of building the transfer matrix between two bases of
% orthogonalized eigenproducts using the analytical formula.
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
MeshFiles = ["Meshes/TOSCA/wolf0.off";
             "Meshes/TOSCA/wolf2.off"];
N = 2;
K = 4;
OrthoThresh = 1e-5;


%% Load the meshes and solve eigenproblem
fprintf("Loading the meshes... ");
tic;
Src = mesh.init(MeshFiles(1));
Trg = mesh.init(MeshFiles(2));
toc;

fprintf("Solve the eigenproblem... ");
tic;
[Src.S, ~, Src.A] = mesh.proc.FEM_higher(Src, 1, 'Dirichlet');
[Src.Phi, Src.Lambda] = eigs(Src.S, Src.A, K, -1e-5);
Src.Lambda = diag(Src.Lambda);

[Trg.S, ~, Trg.A] = mesh.proc.FEM_higher(Trg, 1, 'Dirichlet');
[Trg.Phi, Trg.Lambda] = eigs(Trg.S, Trg.A, K, -1e-5);
Trg.Lambda = diag(Trg.Lambda);
toc;


%% Compute eigenproducts and orthogonalize
fprintf("Computing eigenproducts... ");
tic;
Src.PolyPhi = eigprods(Src, K - 1, N, false);
Trg.PolyPhi = eigprods(Trg, K - 1, N, false);
toc;

fprintf("Orthogonalizing... ");
tic;
[Src.Q, Src.R] = orthogonalize(Src.A, Src.PolyPhi, OrthoThresh, false);
[Trg.Q, Trg.R] = orthogonalize(Trg.A, Trg.PolyPhi, OrthoThresh, false);
toc;


%% Compute GT functional map
C = Trg.Phi' * Trg.A * Src.Phi;
CTilde = Cp_of_C(C, N, size(C, 1) - 1, size(C, 2) - 1, Src.Phi(1, 1));

%% Compute orthotransfer map
O = zeros(size(Trg.Q, 2), size(Src.Q, 2));
O(1:size(C, 1), 1:size(C, 2)) = C;
U = Trg.Q \ (Src.PolyPhi); % N.R*CTilde(:,M.idxs); % N.Q\(N.PolyPhi*CTilde(:,M.idxs)); % 
for idx = K+1:(size(Src.Q, 2))
    O(:, idx) = (1 ./ Src.R(idx,idx)) .* ( U(:, idx) - O(:, 1:idx-1) * Trg.R(1:idx-1, idx));
    O(1:K, idx) = 0;
end
% O = Trg.R * CTilde / Src.R;


%% Plot meshes and O
subplot(2, 2, 1);
utils.plot_scalar_map(Src, zeros(Src.n, 1), utils.cmaps.constant(1), true, 'interp', false);

subplot(2, 2, 2);
utils.plot_scalar_map(Trg, zeros(Trg.n, 1), utils.cmaps.constant(1), true, 'interp', false);
view(-10, 5);

ax = subplot(2, 2, 3);
imagesc(O);
axis image;
colormap(ax, utils.cmaps.bwr);
Lim = max(abs(O), [], 'all');
caxis([-Lim, Lim]);
colorbar;
title("Analytic");


ax = subplot(2, 2, 4);
imagesc(Trg.Q' * Trg.A * Src.Q);
axis image;
colormap(ax, utils.cmaps.bwr);
Lim = max(abs(Trg.Q' * Trg.A * Src.Q), [], 'all');
caxis([-Lim, Lim]);
colorbar;
title("Ground Truth");





























































