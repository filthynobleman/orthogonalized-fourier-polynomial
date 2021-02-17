% Plots the visual transfer and the numerical error for a given function
% between two shapes using the following approaches:
%   - K Laplacian eigenfunctions
%   - N * K Laplacian eigenfunctions
%   - Eigenproducts
%   - Orthogonalized eigenproducts with high threshold
%   - Orthogonalized eigenproducts with low threshold and reiteration.
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


%% Paths
addpath('./mexmesh')
addpath(genpath('./geodesics_matlab/'))
addpath(genpath('./AWFTclean_code_3DV2016/'))
addpath(genpath('ZoomOut2Share/'))
addpath(genpath('./Commute_Fmap_basic_demo/'));
addpath(genpath('./flann/'));


%% Parameters
SrcFile = 'Meshes/man_woman_gorilla/Cman3.mat';
TrgFile = 'Meshes/man_woman_gorilla/Cgorilla2.mat';
KSrc = 30;
KTrg = 40;
N = 2;          % Please, keep this 2 or 3
OrthoThresh = 1e-2;
OrthoThreshReiter = 1e-9;
SVDThresh = (KSrc - 2) * 2e-3;
NumLM = 5;
FontSize = 16;
PlotTitles = true;


%% Load the meshes and compute the needed fields
fprintf("Loading the meshes... ");
tic;

Src = load(SrcFile).shape;
Src.VERT = [Src.X, Src.Y, Src.Z];
Src.n = size(Src.VERT, 1);
Src.m = size(Src.TRIV, 1);
% Src = mesh.transform.rotate(Src, 'x', -90);
Src = mesh.transform.normalize(Src);
LMIdx = datasample(1:length(Src.cor), NumLM, 'Replace', false);
LMIdx = sort(reshape(LMIdx, NumLM, 1));
Src.LM = Src.cor(LMIdx);

Trg = load(TrgFile).shape;
Trg.VERT = [Trg.X, Trg.Y, Trg.Z];
Trg.n = size(Trg.VERT, 1);
Trg.m = size(Trg.TRIV, 1);
% Trg = mesh.transform.rotate(Trg, 'x', -90);
Trg = mesh.transform.normalize(Trg);
Trg.LM = Trg.cor(LMIdx);
toc;

fprintf("Solve the eigenproblem... ");
tic;
[Src.S, ~, Src.A] = mesh.proc.FEM_higher(Src, 1, 'Dirichlet');
[Trg.S, ~, Trg.A] = mesh.proc.FEM_higher(Trg, 1, 'Dirichlet');

[Src.Phi, Src.Lambda] = eigs(Src.S, Src.A, 200, -1e-5);
Src.Lambda = diag(Src.Lambda);

[Trg.Phi, Trg.Lambda] = eigs(Trg.S, Trg.A, 200, -1e-5);
Trg.Lambda = diag(Trg.Lambda);
toc;


%% Initialize a function to transfer and plot it
% F = Src.Y;
% MaxF = max(F);
% MinF = min(F);
% LenF = MaxF - MinF;
% FNorm = (F - MinF) ./ LenF;
% f = sin(4 .* pi .* FNorm);
hk = mesh.proc.heat_kernel(Src, 1e-3, 50);
f = hk(:, 386);
% f = Src.X;
% f = double(Src.Y < 0.66);
Lim = max(abs(f));

subplot(2, 4, 1);
utils.plot_scalar_map(Src, f, utils.cmaps.bwr, true);
if PlotTitles
    title('Source', ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Compute NK functional map and transfer
CNK = build_funmap_largedesc(Src, Trg, N * KSrc, N * KTrg, Src.LM, Trg.LM);
Emb = Src.Phi(:, 1:N*KSrc)' * (Src.A * f);
FTran = Trg.Phi(:, 1:N*KTrg) * (CNK * Emb);
% FTran = (((FTran - min(FTran)) ./ (max(FTran) - min(FTran))) .* LenF) + MinF;
Diff = f(Src.cor) - FTran(Trg.cor);
ErrNK = mean(abs(Diff)) / length(Src.cor);
Lim = max(Lim, max(abs(FTran)));

subplot(2, 4, 3);
utils.plot_scalar_map(Trg, FTran, utils.cmaps.bwr, true);
if PlotTitles
    TitleString = {"$NK$ Eigs";
                   sprintf("Err: %.3f", ErrNK)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Compute K functional map and transfer
C = build_funmap_largedesc(Src, Trg, KSrc, KTrg, Src.LM, Trg.LM);
Emb = Src.Phi(:, 1:KSrc)' * (Src.A * f);
FTran = Trg.Phi(:, 1:KTrg) * (C * Emb);
% FTran = (((FTran - min(FTran)) ./ (max(FTran) - min(FTran))) .* LenF) + MinF;
Diff = f(Src.cor) - FTran(Trg.cor);
ErrK = mean(abs(Diff)) / length(Src.cor);
Lim = max(Lim, max(abs(FTran)));

subplot(2, 4, 4);
utils.plot_scalar_map(Trg, FTran, utils.cmaps.bwr, true);
if PlotTitles
    TitleString = {"$K$ Eigs";
                   sprintf("Err: %.3f", ErrK)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Transfer with eigenproducts
Src.PolyPhi = eigprods(Src, KSrc - 1, N, false);
Trg.PolyPhi = eigprods(Trg, KTrg - 1, N, false);
CTilde = Cp_of_C(C, N, size(C, 1) - 1, size(C, 2) - 1, Src.Phi(1, 1));

[~, S2, V2] = svd(sqrt(Src.A) * Src.PolyPhi);
singular_values = diag(S2);
q_Src = sum(singular_values ./ singular_values(1) > SVDThresh);
tmp = S2 * V2';
projection = tmp(1:q_Src, :); 
clear V2; clear S2; clear tmp;
PinvPoly = pinv(sqrt(Src.A) * Src.PolyPhi);
PinvProj = pinv(projection);
Emb = PinvProj * (projection * (PinvPoly * (sqrt(Src.A) * f)));
FTran = Trg.PolyPhi * (CTilde * Emb);
% FTran = (((FTran - min(FTran)) ./ (max(FTran) - min(FTran))) .* LenF) + MinF;
Diff = f(Src.cor) - FTran(Trg.cor);
ErrPoly = mean(abs(Diff)) / length(Src.cor);
Lim = max(Lim, max(abs(FTran)));

subplot(2, 4, 5);
utils.plot_scalar_map(Trg, FTran, utils.cmaps.bwr, true);
if PlotTitles
    TitleString = {"Poly";
                   sprintf("Err: %.3f", ErrPoly)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Transfer with orthogonal basis
Src.Q = orthogonalize(Src.A, Src.PolyPhi, OrthoThresh, false);
Trg.Q = orthogonalize(Trg.A, Trg.PolyPhi, OrthoThresh, false);

T12 = flann_search(CTilde * Src.PolyPhi', (CTilde * CTilde') * Trg.PolyPhi', 1, struct());
Pi21 = sparse(1:Trg.n, T12, 1, Trg.n, Src.n);
O = Trg.Q' * Trg.A * Pi21 * Src.Q;

Emb = Src.Q' * (Src.A * f);
FTran = Trg.Q * (O * Emb);
% FTran = (((FTran - min(FTran)) ./ (max(FTran) - min(FTran))) .* LenF) + MinF;
Diff = f(Src.cor) - FTran(Trg.cor);
ErrOrtho = mean(abs(Diff)) / length(Src.cor);
Lim = max(Lim, max(abs(FTran)));

subplot(2, 4, 6);
utils.plot_scalar_map(Trg, FTran, utils.cmaps.bwr, true);
if PlotTitles
    TitleString = {"Ortho";
                   sprintf("Err: %.3f", ErrOrtho)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Transfer with stressed orthogonal basis
Src.QS = orthogonalize(Src.A, Src.PolyPhi, OrthoThreshReiter, true);
Trg.QS = orthogonalize(Trg.A, Trg.PolyPhi, OrthoThreshReiter, true);

OS = Trg.QS' * Trg.A * Pi21 * Src.QS;

Emb = Src.QS' * (Src.A * f);
FTran = Trg.QS * (OS * Emb);
% FTran = (((FTran - min(FTran)) ./ (max(FTran) - min(FTran))) .* LenF) + MinF;
Diff = f(Src.cor) - FTran(Trg.cor);
ErrOrthoReiter = mean(abs(Diff)) / length(Src.cor);
Lim = max(Lim, max(abs(FTran)));

subplot(2, 4, 7);
utils.plot_scalar_map(Trg, FTran, utils.cmaps.bwr, true);
if PlotTitles
    TitleString = {"Ortho Reiter";
                   sprintf("Err: %.3f", ErrOrthoReiter)};
    title(TitleString, ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end



%% Caxis all
for i = 1:7
    subplot(2, 4, i);
    caxis([-Lim, Lim]);
    %colorbar;
end


%% Plot the ground truth
subplot(2, 4, 2);
% FGT = zeros(Trg.n, 1);
% FGT(Trg.cor) = f(Src.cor);
utils.plot_scalar_map(Trg, zeros(Trg.n, 1), white, true, 'interp', false);
% utils.plot_scalar_map(Trg, Pi21 * f, utils.cmaps.bwr, true);
% caxis([-Lim, Lim]);
FSample = f(Src.cor);
CM = utils.cmaps.bwr;
FSampleNorm = (FSample + Lim) ./ (2 * Lim);% (max(f) - min(f));
FSampleIdx = min(round(FSampleNorm .* (length(CM) - 1)) + 1, length(CM));
FCols = CM(FSampleIdx, :);
hold on;
TrgPlot = mesh.transform.rotate(Trg, 'z', 30);
scatter3(TrgPlot.X(Trg.cor), TrgPlot.Y(Trg.cor), TrgPlot.Z(Trg.cor), 72, FCols, 'filled');
if PlotTitles
    title('Ground Truth', ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end
hold off;




























































