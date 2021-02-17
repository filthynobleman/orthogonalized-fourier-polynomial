% Transfer the XYZ coordinates function between two textured spheres 
% by using the standard eigenfunctions approach, both with K and N*K 
% eigenfunctions, and compare it to the result obtained using the 
% eigenproducts and the orthogonal basis method. The orthogonal basis, 
% here, is computed using both an high threshold and a low threshold with 
% reiteration.
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
SrcFile = 'Meshes/texturedMeshes/t1_donut.obj';
TrgFile = 'Meshes/texturedMeshes/t0_donut.obj';
KSrc = 27;
KTrg = 30;
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
Trg = mesh.init(TrgFile);
Src.VERT = Src.VERT./sqrt(sum(mesh.proc.tri_areas(Src)));
Trg.VERT = Trg.VERT./sqrt(sum(mesh.proc.tri_areas(Trg)));
p2p_gt_s2ts = knnsearch(Src.VERT,Trg.VERT);
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
f = Src.VERT;
g = f(p2p_gt_s2ts);

%% Plot the original shape
fig = figure;
fig.WindowState = 'maximized';


subplot(2, 4, 1);
utils.plot_scalar_map(Src, zeros(Src.n, 1), white, true, 'interp', false);
title("Source", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize);

  subplot(2, 4, 2);
utils.plot_scalar_map(Trg, zeros(Trg.n, 1), white, true, 'interp', false);
title("Target", ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize);

  
%% Compute NK functional map and transfer
% CNK = build_funmap(Src, Trg, N * KSrc, N * KTrg);
CNK = Trg.Phi(:, 1:N*KTrg)\Src.Phi(p2p_gt_s2ts, 1:N*KSrc);
Emb = Src.Phi(:, 1:N*KSrc)' * (Src.A * f);
FTran = Trg.Phi(:, 1:N*KTrg) * (CNK * Emb);
Diff = g - FTran;
ErrNK = sqrt(Diff' * Trg.A * Diff) / sqrt(g' * Trg.A * g);
New = Trg;
New.VERT = FTran;

subplot(2, 4, 3);
Den = vecnorm(FTran, 2, 2);
Err = vecnorm(FTran - Trg.VERT, 2, 2) ./ Den;
utils.plot_scalar_map(New, Err, flipud(hot(256)), true, 'interp', false);
TitleString = {"$NK$ Eigs"; 
               sprintf("Err: %.3f", sqrt((Err' * Trg.A * Err) / (Den' * Trg.A * Den)))};
title(TitleString, ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize);
colorbar;

%% Compute K functional map and transfer
% C = build_funmap(Src, Trg, KSrc, KTrg);
C = Trg.Phi(:, 1:KTrg)\Src.Phi(p2p_gt_s2ts, 1:KSrc);

Emb = Src.Phi(:, 1:KSrc)' * (Src.A * f);
FTran = Trg.Phi(:, 1:KTrg) * (C * Emb);
Diff = g - FTran;
ErrNK = sqrt(Diff' * Trg.A * Diff) / sqrt(g' * Trg.A * g);
New = Trg;
New.VERT = FTran;

subplot(2, 4, 4);
Den = vecnorm(FTran, 2, 2);
Err = vecnorm(FTran - Trg.VERT, 2, 2) ./ Den;
utils.plot_scalar_map(New, Err, flipud(hot(256)), true, 'interp', false);
TitleString = {"$K$ Eigs"; 
               sprintf("Err: %.3f", sqrt((Err' * Trg.A * Err) / (Den' * Trg.A * Den)))};
title(TitleString, ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize);
colorbar;

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
Diff = g - FTran;
ErrNK = sqrt(Diff' * Trg.A * Diff) / sqrt(g' * Trg.A * g);
New = Trg;
New.VERT = FTran;

subplot(2, 4, 5);
Den = vecnorm(FTran, 2, 2);
Err = vecnorm(FTran - Trg.VERT, 2, 2) ./ Den;
utils.plot_scalar_map(New, Err, flipud(hot(256)), true, 'interp', false);
TitleString = {"Poly"; 
               sprintf("Err: %.3f", sqrt((Err' * Trg.A * Err) / (Den' * Trg.A * Den)))};
title(TitleString, ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize);
colorbar;


%% Transfer with orthogonal basis
Src.Q = orthogonalize(Src.A, Src.PolyPhi, OrthoThresh, false);
Trg.Q = orthogonalize(Trg.A, Trg.PolyPhi, OrthoThresh, false);

T12 = flann_search(CTilde * Src.PolyPhi', (CTilde * CTilde') * Trg.PolyPhi', 1, struct());
Pi21 = sparse(1:Trg.n, T12, 1, Trg.n, Src.n);
O = Trg.Q' * Trg.A * Pi21 * Src.Q;

Emb = Src.Q' * (Src.A * f);
FTran = Trg.Q * (O * Emb);
Diff = g - FTran;
ErrNK = sqrt(Diff' * Trg.A * Diff) / sqrt(g' * Trg.A * g);
New = Trg;
New.VERT = FTran;

subplot(2, 4, 6);
Den = vecnorm(FTran, 2, 2);
Err = vecnorm(FTran - Trg.VERT, 2, 2) ./ Den;
utils.plot_scalar_map(New, Err, flipud(hot(256)), true, 'interp', false);
TitleString = {"Ortho"; 
               sprintf("Err: %.3f", sqrt((Err' * Trg.A * Err) / (Den' * Trg.A * Den)))};
title(TitleString, ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize);
colorbar;


%% Transfer with stressed orthogonal basis
Src.QS = orthogonalize(Src.A, Src.PolyPhi, OrthoThreshReiter, true);
Trg.QS = orthogonalize(Trg.A, Trg.PolyPhi, OrthoThreshReiter, true);

OS = Trg.QS' * Trg.A * Pi21 * Src.QS;

Emb = Src.QS' * (Src.A * f);
FTran = Trg.QS * (OS * Emb);
Diff = g - FTran;
ErrNK = sqrt(Diff' * Trg.A * Diff) / sqrt(g' * Trg.A * g);
New = Trg;
New.VERT = FTran;

subplot(2, 4, 7);
Den = vecnorm(FTran, 2, 2);
Err = vecnorm(FTran - Trg.VERT, 2, 2) ./ Den;
utils.plot_scalar_map(New, Err, flipud(hot(256)), true, 'interp', false);
TitleString = {"Ortho Reiter"; 
               sprintf("Err: %.3f", sqrt((Err' * Trg.A * Err) / (Den' * Trg.A * Den)))};
title(TitleString, ...
      'Interpreter', 'LaTeX', ...
      'FontSize', FontSize);
colorbar;


















































