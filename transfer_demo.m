%% Clear the workspace
clear all
close all
clc


%% Parametrers
KTrg = 39;
KSrc = 29;
SVDThresh = (KSrc - 1) * 2e-3;
OrthoThresh = 1e-9;
WhichO = 1;         % 1: Induce P2P mapping using functional map
                    % 2: Analytical formulation
                    % 3: Ground truth P2P mapping

%% Reproducibility
rng(0);


%% Compute eigenfunctions
fprintf("Loading the meshes... ");
tic;
Trg = mesh.init('Meshes/tr_reg/tr_reg_000.ply');
Src = mesh.init('Meshes/tr_reg/tr_reg_020.ply');
% Trg = mesh.init('Meshes/tr_reg/tr_reg_094.ply');
% Src = mesh.init('Meshes/tr_reg/tr_reg_019.ply');
toc;

fprintf("Computing laplacian eigenfunctions... ");
tic;
[Trg.S, ~, Trg.A] = mesh.proc.FEM_higher(Trg, 1, 'Dirichlet');
[Trg.Phi, Trg.Lambda] = eigs(Trg.S, Trg.A, KTrg + 1, -1e-5);
Trg.Lambda = diag(Trg.Lambda);

[Src.S, ~, Src.A] = mesh.proc.FEM_higher(Src, 1, 'Dirichlet');
[Src.Phi, Src.Lambda] = eigs(Src.S, Src.A, KSrc + 1, -1e-5);
Src.Lambda = diag(Src.Lambda);
toc;

fprintf("Computing functional map... ");
tic;
C = Trg.Phi' * (Trg.A * Src.Phi);
toc;


%% Compute eigenproducts
fprintf("Computing eigenproducts... ");
tic;
Trg.PolyPhi = eigprods(Trg, KTrg, 2, false);
Src.PolyPhi = eigprods(Src, KSrc, 2, false);
toc;

fprintf("Extending functional map... ");
tic;
CTilde = eigprod_funmap(C, Trg.Phi(:, 1));
toc;


%% Compute orthogonal
fprintf("Eigenproducts orthonormalization... ");
tic;
[Trg.Q, Trg.R] = orthogonalize(Trg.A, Trg.PolyPhi, OrthoThresh, true);
[Src.Q, Src.R] = orthogonalize(Src.A, Src.PolyPhi, OrthoThresh, true);
toc;

fprintf("Computing orthogonal extended map... ");
tic;
if WhichO == 1
    T12 = knnsearch(Src.PolyPhi * CTilde', Trg.PolyPhi * (CTilde * CTilde'));
    Pi21 = sparse(1:Trg.n, T12, 1, Trg.n, Src.n);
    O = Trg.Q' * Trg.A * Pi21 * Src.Q;
elseif WhichO == 2
    O = Trg.R * CTilde / Src.R;
elseif WhichO == 3
    O = Trg.Q' * Trg.A * speye(Trg.n, Src.n) * Src.Q;
else
    error("Invalid value (%d) for parameter WhichO", WhichO);
end
toc;



%% Compare transfers
f = Src.Z;
% f = mesh.proc.heat_kernel(M2, 1e-3, 200);
% f = f(:, 1);
% f = M2.VERT;
% f = 1 + cos(35 .* (M2.Phi(:, 2) - min(M2.Phi(:, 2))) ./ (max(M2.Phi(:, 2)) - min(M2.Phi(:, 2))));

fig = figure;
fig.WindowState = 'maximized';

subplot(2, 6, [1 2 3]);
utils.plot_scalar_map(Src, f, utils.cmaps.bwr);
colorbar;
title("Source", ...
      'FontSize', 24);

subplot(2, 6, [4 5 6]);
utils.plot_scalar_map(Trg, f, utils.cmaps.bwr);
colorbar;
title("Ground Truth", ...
      'FontSize', 24);




%% Transfer using eigenfunctions
fenc = Src.Phi' * Src.A * f;
fenctran = C * fenc;
ftran = Trg.Phi * fenctran;
d = f - ftran;
e = sqrt(d' * Trg.A * d) / (f' * Trg.A * f);

subplot(2, 6, [7 8]);
utils.plot_scalar_map(Trg, ftran(1:Trg.n), utils.cmaps.bwr);
colorbar
title(sprintf("Eigenfunctions Err: %.3f", e), ...
      'FontSize', 24);




%% Transfer using eigenproducts
[~, S2, V2] = svd(sqrt(Src.A) * Src.PolyPhi);
singular_values = diag(S2);
q_Src = sum(singular_values ./ singular_values(1) > SVDThresh);
tmp = S2 * V2';
projection = tmp(1:q_Src, :); 
clear V2; clear S2; clear tmp;

fenc = pinv(projection) * (projection * ((pinv(Src.PolyPhi)) * f));
fenctran = CTilde * fenc;
ftran = Trg.PolyPhi * fenctran;
d = f - ftran;
e = sqrt(d' * Trg.A * d) / (f' * Trg.A * f);

subplot(2, 6, [9 10]);
utils.plot_scalar_map(Trg, ftran(1:Trg.n), utils.cmaps.bwr);
colorbar
title(sprintf("Eigenproducts Err: %.3f", e), ...
      'FontSize', 24);


%% Transfer using orthogonal
fenc = Src.Q' * Src.A * f;
fenctran = O * fenc;
ftran = Trg.Q * fenctran;
d = f - ftran;
e = sqrt(d' * Trg.A * d) / (f' * Trg.A * f);

subplot(2, 6, [11 12]);
utils.plot_scalar_map(Trg, ftran(1:Trg.n), utils.cmaps.bwr);
colorbar
title(sprintf("Orthogonal Err: %.3f", e), ...
      'FontSize', 24);













































