% Example showing the error induced by the analytic expression for the
% orthogonal basis transfer map.
%
% This script reproduces Fig. 10 (right) of the reference paper.
%
% This code is part of the implementation of the paper "Orthogonalized 
% Fourier Polynomials for Signal Approximation and Transfer", authored by
% Filippo Maggioli, Simone Melzi, Michael Bronstein, Maks Ovsjanikov and
% Emanuele Rodolà.

%% Clear workspace
clear all
close all
clc


%% Parameters
K = 5:1:20;
S2TRatio = [0.5, 1, 2];
OrthoThresh = 1e-6;
MeshFiles = {'Meshes/tr_reg/tr_reg_000.ply', 'Meshes/tr_reg/tr_reg_019.ply'};
FontSize = 24;



%% Load the mesh and compute eigenfunctions
fprintf("Loading the meshes and solving eigenproblem... ");
tic;
Src = mesh.init(MeshFiles{1});
[Src.S, ~, Src.A] = mesh.proc.FEM_higher(Src, 1, 'Dirichlet');
[Src.Phi, Src.Lambda] = eigs(Src.S, Src.A, max(K), -1e-5);
Src.Lambda = diag(Src.Lambda);

Trg = mesh.init(MeshFiles{2});
[Trg.S, ~, Trg.A] = mesh.proc.FEM_higher(Trg, 1, 'Dirichlet');
[Trg.Phi, Trg.Lambda] = eigs(Trg.S, Trg.A, ceil(max(K) * max(S2TRatio)), -1e-5);
Trg.Lambda = diag(Trg.Lambda);
toc;


%% Compute error
ErrReal =   zeros(length(K), length(S2TRatio));

for i = 1:length(K)
    for j = 1:length(S2TRatio)
        fprintf("Iteration %d out of %d... ", (i - 1) * length(S2TRatio) + j, length(K) * length(S2TRatio));
        tic;

        KSrc = K(i);
        KTrg = ceil(KSrc * S2TRatio(j));

        % Build Functional map and extension
        C = Trg.Phi(:, 1:KTrg)'* (Trg.A * Src.Phi(:, 1:KSrc));
        CTilde = eigprod_funmap(C, Src.Phi(1, 1));

        % Compute eigenproducts
        Src.PolyPhi = eigprods(Src, KSrc - 1, 2, false);
        Trg.PolyPhi = eigprods(Trg, KTrg - 1, 2, false);

        % Orthogonalize
        [Src.Q, Src.R] = orthogonalize(Src.A, Src.PolyPhi, OrthoThresh, true);
        [Trg.Q, Trg.R] = orthogonalize(Trg.A, Trg.PolyPhi, OrthoThresh, true);

        % Compute O
        O = Trg.R * (CTilde / Src.R);

        % Compute errors
        ErrReal(i, j) = norm(Src.Q - Trg.Q * O, 'fro');


        toc;
    end
end


%% Plot
fig = figure;
fig.WindowState = 'maximized';

hold on;
plot(K, ErrReal(:, 1), '-o', ...
     'LineWidth', 3);
plot(K, ErrReal(:, 2), '-^', ...
     'LineWidth', 3);
plot(K, ErrReal(:, 3), '-*', ...
     'LineWidth', 3);
grid; grid minor;
set(gca, 'YScale', 'log');

xlabel('Number of Eigenfunctions', ...
       'Interpreter', 'LaTeX', ...
       'FontSize', FontSize);
ylabel('Error', ...
       'Interpreter', 'LaTeX', ...
       'FontSize', FontSize);
Labels = sprintf('$|\\Phi| / |\\Psi| = %.1f$&', S2TRatio);
Labels = split(Labels, '&');
legend(Labels(1:end-1), ...
       'Interpreter', 'LaTeX', ...
       'FontSize', FontSize, ...
       'Location', 'northwest');
































































