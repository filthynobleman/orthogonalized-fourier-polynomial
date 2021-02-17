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


%% Parameters
SrcFile = 'Meshes/SHREC/32.off';
TrgFile = 'Meshes/SHREC/10.off';
KSrc = 30;
KTrg = 20;
N = 2;          % Please, keep this 2 or 3
OrthoThresh = 1e-2;
OrthoThreshReiter = 1e-9;
SVDThresh = (KSrc - 2) * 2e-3;
FontSize = 24;
PlotTitles = true;


%% Load the meshes and compute the needed fields
fprintf("Loading the meshes... ");
tic;
Src = mesh.init(SrcFile);
Src = mesh.transform.rotate(Src, 'z', -45);
Src = mesh.transform.normalize(Src);
Trg = mesh.init(TrgFile);
Trg = mesh.transform.rotate(Trg, 'x', 90);
Trg = mesh.transform.normalize(Trg);
toc;

fprintf("Loading GT correspondence... ");
tic;
Src.GT = load('Meshes/SHREC/10_32.mat').gt_matches;
Pi21 = sparse(1:Trg.n, Src.GT, 1, Trg.n, Src.n);
toc;

fprintf("Solve the eigenproblem... ");
tic;
[Src.S, ~, Src.A] = mesh.proc.FEM_higher(Src, 1, 'Dirichlet');
[Trg.S, ~, Trg.A] = mesh.proc.FEM_higher(Trg, 1, 'Dirichlet');

[Src.Phi, Src.Lambda] = eigs(Src.S, Src.A, N * KSrc, -1e-5);
Src.Lambda = diag(Src.Lambda);

[Trg.Phi, Trg.Lambda] = eigs(Trg.S, Trg.A, N * KTrg, -1e-5);
Trg.Lambda = diag(Trg.Lambda);
toc;


%% Initialize a function to transfer and plot it
F = Src.Phi(:, 8);
NormF = (F - min(F)) ./ (max(F) - min(F));
f = sin(2 * pi * NormF);
% f = sin(2 .* pi .* Src.X);
% hk = mesh.proc.heat_kernel(Src, 1e-3, 200);
% f = hk(:, 6952) + hk(:, 3651);
% f = Src.X .* Src.Y;
MinF = min(f);
MaxF = max(f);
LenF = MaxF - MinF;
Lim = max(abs(f));

subplot(2, 4, 1);
utils.plot_scalar_map(Src, f, utils.cmaps.bwr, true);
if PlotTitles
    title('Source', ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end
subplot(2, 4, 2);
utils.plot_scalar_map(Trg, f(Src.GT), utils.cmaps.bwr, true);
if PlotTitles
    title('Ground Truth', ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Compute NK functional map and transfer
% CNK = build_funmap(Src, Trg, N * KSrc, N * KTrg);
CNK = Trg.Phi(:, 1:N*KTrg)' * Trg.A * Pi21 * Src.Phi(:, 1:N*KSrc);
Emb = Src.Phi(:, 1:N*KSrc)' * (Src.A * f);
FTran = Trg.Phi(:, 1:N*KTrg) * (CNK * Emb);
Diff = f(Src.GT) - FTran;
ErrNK = sqrt(Diff' * Trg.A * Diff) / sqrt(f(Src.GT)' * Trg.A * f(Src.GT));
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
% C = build_funmap(Src, Trg, KSrc, KTrg);
C = Trg.Phi(:, 1:KTrg)' * Trg.A * Pi21 * Src.Phi(:, 1:KSrc);
Emb = Src.Phi(:, 1:KSrc)' * (Src.A * f);
FTran = Trg.Phi(:, 1:KTrg) * (C * Emb);
% FTran = (((FTran - min(FTran)) ./ (max(FTran) - min(FTran))) .* LenF) + MinF;
Diff = f(Src.GT) - FTran;
ErrK = sqrt(Diff' * Trg.A * Diff) / sqrt(f(Src.GT)' * Trg.A * f(Src.GT));
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
Diff = f(Src.GT) - FTran;
ErrPoly = sqrt(Diff' * Trg.A * Diff) / sqrt(f(Src.GT)' * Trg.A * f(Src.GT));
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
Pi21Comp = sparse(1:Trg.n, T12, 1, Trg.n, Src.n);
O = Trg.Q' * Trg.A * Pi21Comp * Src.Q;

Emb = Src.Q' * (Src.A * f);
FTran = Trg.Q * (O * Emb);
% FTran = (((FTran - min(FTran)) ./ (max(FTran) - min(FTran))) .* LenF) + MinF;
Diff = f(Src.GT) - FTran;
ErrOrtho = sqrt(Diff' * Trg.A * Diff) / sqrt(f(Src.GT)' * Trg.A * f(Src.GT));
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
Diff = f(Src.GT) - FTran;
ErrOrthoReiter = sqrt(Diff' * Trg.A * Diff) / sqrt(f(Src.GT)' * Trg.A * f(Src.GT));
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
end





























































