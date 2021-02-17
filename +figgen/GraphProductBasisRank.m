% Plot the number of linearly independent eigenproducts against the total
% number of eigenproducts, at varying of the number of eigenfunctions
% involved and the order of the products.
%
% This script reproduces Fig. 6 of the reference paper.
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
MeshFiles = ["Meshes/tr_reg/tr_reg_000.ply";
             "Meshes/bunny.off";
             "Meshes/cat0.off";
             "Meshes/texturedMeshes/t1_donut.obj"];
N = 1:5;
K = [   1,     2,      3,      4,      5;
       10,    20,     30,     40,     50;
        1,     5,     10,     15,     20;
        1,     5,     10,     15,     20;
        1,     3,      6,      9,     12];
FontSize = 24;


%% Reproducibility
rng(0);


%% Load the mesh and compute stiff/mass
fprintf("Loading the meshes... ");
tic;
M = cell(length(MeshFiles), 1);
for i = 1:length(MeshFiles)
    M{i} = mesh.init(MeshFiles(i));
end
toc;

fprintf("Computing laplacian... ");
tic;
for i = 1:length(M)
    [M{i}.S, M{i}.A] = mesh.proc.FEM_higher(M{i}, 1, 'Dirichlet');
end
toc;

fprintf("Computing max number of eigenfunctions");
tic;
for i = 1:length(M)
    [M{i}.Phi, M{i}.Lambda] = eigs(M{i}.S, M{i}.A, max(K, [], 'all'), -1e-5);
    M{i}.Lambda = diag(M{i}.Lambda);
end
toc


%% Compute basis rank
Ranks = zeros(length(N), size(K, 2), length(M));
for h = 1:length(M)
    for i = 1:length(N)
        for j = 1:size(K, 2)
            n = N(i);
            k = K(i, j);
            fprintf("N = %d, K = %d\n", n, k);
            PolyPhi = eigprods(M{h}, k, n, true);
            Ranks(i, j, h) = rank(PolyPhi);
            clear PolyPhi;
        end
    end
end


%% Save results
save Results.ProdRank.Mean.mat


%% Plot each N
fig3 = figure;
fig3.WindowState = 'maximized';

Ranks = mean(Ranks, 3);

for i = 1:4
    subplot(2, 2, i);
    L = size(K, 2);
    Ref = zeros(L, 1);
    for j = 1:L
        Ref(j) = nchoosek(N(i + 1) + K(i + 1, j), N(i + 1));
    end
    hold on;
    plot(K(i + 1, :), Ranks(i + 1, 1:L), ...
         '-o', ...
         'LineWidth', 3);
    scatter(K(i + 1, :), Ref, ...
         '^', ...
         'LineWidth', 3);
    xlabel("Number of Eigenfunctions", ...
           'Interpreter', 'LaTeX', ...
           'FontSize', FontSize);
    ylabel("Products Basis Rank", ...
           'Interpreter', 'LaTeX', ...
           'FontSize', FontSize);
    legend(sprintf("Products of Order %d", N(i + 1)), sprintf("${{K + %d} \\choose %d}$", N(i + 1), N(i + 1)), ...
           'Interpreter', 'LaTeX', ...
           'FontSize', FontSize, ...
           'Location', 'northwest');
    grid; grid minor;
end
































































