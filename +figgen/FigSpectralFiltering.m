% Example showing the stability of the orthogonalized basis (with respect
% to the eigenproducts) on spectral filtering of the XYZ coordinates
% functions.
%
% This script reproduces Fig. 14 of the reference paper.
%
% This code is part of the implementation of the paper "Orthogonalized 
% Fourier Polynomials for Signal Approximation and Transfer", authored by
% Filippo Maggioli, Simone Melzi, Michael Bronstein, Maks Ovsjanikov and
% Emanuele Rodolà.

%% Clean the workspace
clear all
close all
clc


%% Reproducibility
Seed = 0;
rng(Seed);


%% Parameters
SrcFile = 'Meshes/sheep.off';
KSrc = 15;
N = 3;          % Please, keep this 2 or 3
OrthoThresh = 1e-2;
OrthoThreshReiter = 1e-9;
SVDThresh = (KSrc - 2) * 2e-3;
FontSize = 14;
PlotTitles = true;
tau = 1e-2;

%% Load the meshes and compute the needed fields
fprintf("Loading the meshes... ");
tic;
Src = mesh.init(SrcFile);
Src.VERT = Src.VERT./sqrt(sum(mesh.proc.tri_areas(Src)));
toc;

fprintf("Solve the eigenproblem... ");
tic;
[Src.S, ~, Src.A] = mesh.proc.FEM_higher(Src, 1, 'Dirichlet');
[Src.Phi, Src.Lambda] = eigs(Src.S, Src.A, N * KSrc, -1e-5);
Src.Lambda = diag(Src.Lambda);
toc;


%% Initialize a function to transfer and plot it
f = Src.VERT; % sin(16 .* pi .* Src.X .* Src.Y);

%% Transfer with eigenproducts
Src.PolyPhi = eigprods(Src, KSrc - 1, N, true);
freqs = proc.dirichlet(Src, Src.PolyPhi);
[freqs,indices] = sort(freqs,'ascend'); clear tmp;
basis_poly = Src.PolyPhi;
basis_poly = basis_poly(:,indices); clear indices;
% filter = exp(-0.01.*(abs(-freqs(1)+freqs)))';
filter = exp(-tau*(abs(0:size(basis_poly,2)-1)))';
filter = filter./max(filter);

% coeff = pinv(basis_poly)*f;
coeff = pinv(sqrt(Src.A) * basis_poly) * sqrt(Src.A) * f;
prod = basis_poly*(filter.*coeff);

%% Transfer with stressed orthogonal basis
Src.QS = orthogonalize(Src.A, Src.PolyPhi, OrthoThreshReiter, true);
freqs = proc.dirichlet(Src, Src.QS);
[freqs,indices] = sort(freqs,'ascend'); clear tmp;
basis_our = Src.QS;
basis_our = basis_our(:,indices); clear indices;
% filter = exp(-0.01.*(abs(-freqs(1)+freqs)))';
% filter = exp(-tau*(abs(0:size(basis_our,2)-1)))';
% filter = filter./max(filter);
% figure; plot(filter)

coeff = basis_our' * (Src.A * f);
ortho_forced = basis_our*(filter.*coeff);

%%
New = Src;
figure
subplot(2, 4, 1);
PlotMesh = mesh.transform.rotate(Src, 'x', 170);
PlotMesh = mesh.transform.rotate(PlotMesh, 'y', 25);
utils.plot_scalar_map(PlotMesh, zeros(Src.n, 1), white, true, 'interp', false);
title("Source", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize); %view([0 90])
New.VERT = prod;
subplot(2, 4, 2);
PlotMesh = mesh.transform.rotate(New, 'x', 170);
PlotMesh = mesh.transform.rotate(PlotMesh, 'y', 25);
utils.plot_scalar_map(PlotMesh, zeros(Src.n, 1), white, true, 'interp', false);
title("Poly", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize); %view([0 90])
New.VERT = ortho_forced;
subplot(2, 4, 3);
PlotMesh = mesh.transform.rotate(New, 'x', 170);
PlotMesh = mesh.transform.rotate(PlotMesh, 'y', 25);
utils.plot_scalar_map(PlotMesh, zeros(Src.n, 1), white, true, 'interp', false);
title("Otho", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize); %view([0 90])
subplot(2, 4, 4)
plot(filter,'LineWidth',2)
Filter1 = filter;

%% Transfer with eigenproducts
filter = exp(-0.0025*((abs(size(basis_poly,2)-1)- [0:size(basis_poly,2)-1])))';
% filter = ones(size(basis_poly,2),1)+rand(size(basis_poly,2),1)*0.5; % 
% filter = exp(-tau*(abs(0:size(basis,2)-1)))';
% filter = filter./max(filter);

% coeff = pinv(basis)*f;
coeff = pinv(sqrt(Src.A) * basis_poly) * sqrt(Src.A) * f;
prod = basis_poly*(filter.*coeff);

%% Transfer with stressed orthogonal basis

coeff = basis_our' * (Src.A * f);
ortho_forced = basis_our*(filter.*coeff);

%%
New = Src;
New.VERT = prod;
subplot(2, 4, 5);
PlotMesh = mesh.transform.rotate(Src, 'x', 170);
PlotMesh = mesh.transform.rotate(PlotMesh, 'y', 25);
utils.plot_scalar_map(PlotMesh, zeros(Src.n, 1), white, true, 'interp', false);
title("Source", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize); %view([0 90])
subplot(2, 4, 6);
PlotMesh = mesh.transform.rotate(New, 'x', 170);
PlotMesh = mesh.transform.rotate(PlotMesh, 'y', 25);
utils.plot_scalar_map(PlotMesh, zeros(Src.n, 1), white, true, 'interp', false);
title("Poly", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize); %view([0 90])
New.VERT = ortho_forced;
subplot(2, 4, 7);
PlotMesh = mesh.transform.rotate(New, 'x', 170);
PlotMesh = mesh.transform.rotate(PlotMesh, 'y', 25);
utils.plot_scalar_map(PlotMesh, zeros(Src.n, 1), white, true, 'interp', false);
title("Otho", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize); %view([0 90])
subplot(2, 4, 8)
plot(filter,'LineWidth',2)
Filter2 = filter;


figure;
subplot(2, 1, 1);
plot(Filter1, 'LineWidth', 3);
subplot(2, 1, 2);
plot(Filter2, 'LineWidth', 3);







































