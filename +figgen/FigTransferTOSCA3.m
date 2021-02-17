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
SrcFile = 'Meshes/TOSCA/wolf0.off';
TrgFile = 'Meshes/TOSCA/wolf2.off';
KSrc = 15;
KTrg = 18;
N = 3;          % Please, keep this 2 or 3
OrthoThresh = 1e-2;
OrthoThreshReiter = 1e-9;
SVDThresh = (KSrc - 2) * 2e-3;
FontSize = 24;
PlotTitles = true;


%% Load the meshes and compute the needed fields
fprintf("Loading the meshes... ");
tic;
Src = mesh.init(SrcFile);
% Src = mesh.transform.rotate(Src, 'x', -90);
% Src = mesh.transform.rotate(Src, 'y', 20);
Src = mesh.transform.normalize(Src);
Src.LM = euclidean_fps(Src, 5);
Trg = mesh.init(TrgFile);
% Trg = mesh.transform.rotate(Trg, 'x', -90);
% Trg = mesh.transform.rotate(Trg, 'y', 40);
Trg = mesh.transform.normalize(Trg);
Trg.LM = Src.LM;
toc;

fprintf("Solve the eigenproblem... ");
tic;
[Src.S, ~, Src.A] = mesh.proc.FEM_higher(Src, 1, 'Dirichlet');
[Trg.S, ~, Trg.A] = mesh.proc.FEM_higher(Trg, 1, 'Dirichlet');

[Src.Phi, Src.Lambda] = eigs(Src.S, Src.A, N * KSrc + 5, -1e-5);
Src.Lambda = diag(Src.Lambda);

[Trg.Phi, Trg.Lambda] = eigs(Trg.S, Trg.A, N * KTrg + 5, -1e-5);
Trg.Lambda = diag(Trg.Lambda);
toc;


%% Initialize a function to transfer and plot it
% f = sin(8 .* pi .* Src.X .* Src.Y);
% hk = mesh.proc.heat_kernel(Src, 1e-3, 200);
% f = atan(0.0625 * (diag(hk) ./ max(diag(hk))) ./ (pi / 2));
% f = hk(:, 6952) + hk(:, 3651);
% f = exp(200 * -(Src.X - mean(Src.X)).^2);
% f = sum(mesh.proc.hks(Src, logspace(-8, -3, 6), 200), 2);
% NormX = (Src.X - min(Src.X)) ./ (max(Src.X) - min(Src.X));
% NormY = (Src.Y - min(Src.Y)) ./ (max(Src.Y) - min(Src.Y));
% NormZ = (Src.Z - min(Src.Z)) ./ (max(Src.Z) - min(Src.Z));
% f = sin(2 * pi * NormZ.^2) .* sin(2 * pi * NormY.^6);

% [Funcs, fidxs] = funcset(Src, 20, N * KSrc, 10);
% f = Funcs(:, 45);
f = sum(Src.Phi(:, N * KSrc + 1), 2);

% f = dijkstra_to_all(Src, 1);
% f = round(f .* 10) ./ 10;

% AWFTparam.freqs0 = [1 1];
% AWFTparam.taus = [0.002 0.05]; 
% AWFTparam.angles = [-1 90 180]; % -1 is for the isotropic case
% AWFTparam.curv_smooth = 10;
% AWFTparam.n_eigen = 100;
% AWFTparam.alpha = [ 100];
% AWFTfunctions{1} = 'ShapeIndex';
% AWFTfunctions{2} = 'Fiedler';
% AWFT = compute_AWFTdesc(Src, AWFTfunctions, AWFTparam);
% add =  AWFT;
% add = add./max(abs(add));
% f = add(:, 1);
% WKS = waveKernelSignature_exact(Src.Phi, Src.Lambda, 100);
% WKS_idxs = randperm(100);
% add =  WKS(:, WKS_idxs(1:10));
% add = add ./ max(abs(add));
% f = add;

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
utils.plot_scalar_map(Trg, f, utils.cmaps.bwr, true);
if PlotTitles
    title('Ground Truth', ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
end


%% Compute NK functional map and transfer
CNK = build_funmap(Src, Trg, N * KSrc, N * KTrg, Src.LM, Trg.LM);
Emb = Src.Phi(:, 1:N*KSrc)' * (Src.A * f);
FTran = Trg.Phi(:, 1:N*KTrg) * (CNK * Emb);
% FTran = (((FTran - min(FTran)) ./ (max(FTran) - min(FTran))) .* LenF) + MinF;
Diff = f - FTran;
ErrNK = sqrt(Diff' * Trg.A * Diff) / sqrt(f' * Trg.A * f);
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
C = build_funmap(Src, Trg, KSrc, KTrg, Src.LM, Trg.LM);
Emb = Src.Phi(:, 1:KSrc)' * (Src.A * f);
FTran = Trg.Phi(:, 1:KTrg) * (C * Emb);
% FTran = (((FTran - min(FTran)) ./ (max(FTran) - min(FTran))) .* LenF) + MinF;
Diff = f - FTran;
ErrK = sqrt(Diff' * Trg.A * Diff) / sqrt(f' * Trg.A * f);
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
Diff = f - FTran;
ErrPoly = sqrt(Diff' * Trg.A * Diff) / sqrt(f' * Trg.A * f);
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
Diff = f - FTran;
ErrOrtho = sqrt(Diff' * Trg.A * Diff) / sqrt(f' * Trg.A * f);
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
Diff = f - FTran;
ErrOrthoReiter = sqrt(Diff' * Trg.A * Diff) / sqrt(f' * Trg.A * f);
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





























































