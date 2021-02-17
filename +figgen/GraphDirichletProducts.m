% Computes and plot the error between the first N*K Laplacian eigenvalues
% and the closest Dirichlet energy between the eigenproducts of order K
% involving the first K eigenfunctions. The results are averaged on 4
% shapes.
%
% This script reproduces Fig. 4 of the reference paper.
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
K = 30;
N = 3;
PlotSumHalf = false;
FontSize = 24;


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
    [M{i}.Phi, M{i}.Lambda] = eigs(M{i}.S, M{i}.A, N * K + 1, -1e-5);
    M{i}.Lambda = diag(M{i}.Lambda);
end
toc


%% For each order, store the Dirichlet energies of the products
for j = 1:length(M)
    M{j}.Energies = cell(N, 1);
    for i = 1:N
        fprintf("Dirichlet energies products order %d... ", i);
        tic;
        M{j}.PolyPhi = eigprods(M{j}, K, i, true);
        M{j}.Energies{i} = proc.dirichlet(M{j}, M{j}.PolyPhi)';
        toc;
    end
    M{j} = rmfield(M{j}, 'PolyPhi');
end


%% For each order and each eigenvalue up to N*K compute diff with nearest Dirichlet energy
fprintf("Compute differences... ");
tic;
for j = 1:length(M)
    M{j}.Differences = cell(N, 1);
    for i = 1:N
        M{j}.Differences{i} = min(abs(M{j}.Lambda(1:i*K) - M{j}.Energies{i}'), [], 2);
        M{j}.Differences{i}(2:end) = M{j}.Differences{i}(2:end) ./ M{j}.Lambda(2:i*K);
    end
end
toc;


% %% Compute the sums of eigenvalues
% M.EigSum = cell(N, 1);
% M.EigSum{1} = M.Lambda';
% tmpsum = M.Lambda(2:end)';
% for i = 2:N
%     tmp1 = repmat(M.Lambda(2:end)', 1, K^(i - 1));
%     tmp2 = reshape(repmat(tmpsum, K, 1), 1, K^i);
%     tmpsum = tmp1 + tmp2;
%     M.EigSum{i} = [M.EigSum{i - 1}, tmpsum];
% end
% for i = 2:N
%     idx = sub_index(i, K);
%     M.EigSum{i} = M.EigSum{i}(idx);
% end
% clear tmp1 tmp2 tmpsum idx
% 
% 
% %% For each order, plot energies, eigenvalues sum and eigenvalues sum / 2
% fig = figure;
% fig.WindowState = 'maximized';
% Energies = cell(N, 1);
% for i = 1:N
%     subplot(N, 1, i);
%     hold on;
%     X = 1:length(M.Energies{i});
%     [Energies{i}, idx] = sort(M.Energies{i});
%     p1 = scatter(X, Energies{i}, 'filled');
%     p2 = scatter(X, M.EigSum{i}(idx), 'filled');
%     if PlotSumHalf
%         p3 = scatter(X, M.EigSum{i}(idx) ./ 2, 'filled');
%         uistack(p3, 'bottom');
%     end
%     grid; grid minor;
%     
%     Labels = {'$\mathcal{E}(Q_I)$';
%               '$\sum_{i \in I} \lambda_i$';
%               '$\frac{1}{2}\sum_{i \in I} \lambda_i$'};
%     if ~PlotSumHalf
%         Labels(3) = [];
%     end
%     legend(Labels, ...
%            'Interpreter', 'LaTeX', ...
%            'FontSize', FontSize, ...
%            'Location', 'northwest');
%     title(sprintf("%d Eigenvalues, Order %d", K, i), ...
%           'Interpreter', 'LaTeX', ...
%           'FontSize', FontSize);
%       
%     
%     uistack(p1, 'top');
%     if PlotSumHalf
%         uistack(p3, 'bottom');
%     end
% end


%% For each order, plot the difference between eigenvalues and energies
fig2 = figure;
fig2.WindowState = 'maximized';

for i = 1:N
    Diff = M{1}.Differences{i};
    for j = 2:length(M)
        Diff = Diff + M{j}.Differences{i};
    end
    Diff = Diff ./ length(M);
    
    subplot(N, 1, i);
    X = 1:i*K;
    plot(X, Diff, ...
         'LineWidth', 3);
    grid; grid minor;
    
    xlabel("Eigenvalue Index", ...
           'Interpreter', 'LaTeX', ...
           'FontSize', FontSize);
    ylabel("Relative Distance", ...
           'Interpreter', 'LaTeX', ...
           'FontSize', FontSize);
    title(sprintf("%d Eigenvalues, Order %d", K, i), ...
          'Interpreter', 'LaTeX', ...
          'FontSize', FontSize);
    ylim([0, 0.1]);
end































































