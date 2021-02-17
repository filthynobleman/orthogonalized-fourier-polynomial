% Example showing the stability of the orthogonalized basis (with respect
% to the eigenproducts) on spectral filtering of the XYZ coordinates
% functions.
%
% This script reproduces Fig. 15 of the reference paper.
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

%% Load the meshes and compute the needed fields
fprintf("Loading the meshes... ");
tic;
Src = mesh.init(SrcFile);
% Src = mesh.transform.rotate(Src, 'x', 90);
% Src = mesh.transform.rotate(Src, 'z', -22.5);
% Src = mesh.transform.rotate(Src, 'y', -10);
Src.VERT = Src.VERT./sqrt(sum(mesh.proc.tri_areas(Src)));
% Src = mesh.transform.normalize(Src);
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

% coeff = pinv(basis_poly)*f;
coeff_poly = pinv(sqrt(Src.A) * basis_poly) * sqrt(Src.A) * f;

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

coeff_ours = basis_our' * (Src.A * f);

%%
Prod = Src;
Ours = Src;
w = 0.35;
idx_poly = 4; %8; % randi(size(basis_poly,2),1);
idx_ours = 2; %randi(size(basis_our,2),1);
% [~,idx_poly] = max(sum(coeff_poly,2));
% [~,idx_ours] = max(sum(coeff_ours,2));
def_poly = zeros(length(coeff_poly),1);
def_poly(idx_poly) = w.*coeff_poly(idx_poly);
def_ours = zeros(length(coeff_ours),1);
def_ours(idx_ours) = w.*coeff_ours(idx_ours);
nstep = 6;
fig = figure;
fig.WindowState = 'maximized';
for i = 1:nstep
    
    prod_ = basis_poly*(coeff_poly+(((i-1)/(nstep-1)).*def_poly));
    ours_ = basis_our*(coeff_ours+(((i-1)/(nstep-1)).*def_ours));

    Prod.VERT = prod_;
    Ours.VERT = ours_;

    subplot(2, nstep, i);
    ProdPlot = mesh.transform.rotate(Prod, 'x', 170);
    ProdPlot = mesh.transform.rotate(ProdPlot, 'y', 25);
    utils.plot_scalar_map(ProdPlot, zeros(Src.n, 1), white, true, 'interp', false);
    title("Poly", ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize); %view(0, 90);
   
    subplot(2, nstep, nstep+i);
    ProdPlot = mesh.transform.rotate(Ours, 'x', 170);
    ProdPlot = mesh.transform.rotate(ProdPlot, 'y', 25);
    utils.plot_scalar_map(ProdPlot, zeros(Src.n, 1), white, true, 'interp', false);
    title("Ortho", ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize); %view(0, 90);
end
































