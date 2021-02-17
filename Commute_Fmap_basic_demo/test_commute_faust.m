% This code implements a basic version of the algorithm described in:
% 
% Informative Descriptor Preservation via Commutativity for Shape Matching,
% Dorian Nogneng and Maks Ovsjanikov, Proc. Eurographics 2017
% 
% To try it, simply run this file in MATLAB. This should produce
% a map (correspondence) between a pair of meshes from the FAUST dataset,
% and create an image that visualizes this correspondence.
% 
% This code was written by Etienne Corman and modified by Maks Ovsjanikov.

clear all; close all; 

addpath(genpath(pwd));

%% Load meshes and compute Laplacian eigendecomposition

% Number of basis vectors for computing the functional map.
% Larger is usually better (more accurate results) but somewhat slower.
numEigsSrc = 60;
numEigsTar = 60;

%%
[X, T] = readOff('./Mesh/tr_reg_089');
fprintf('reading the source shape...');tic;
Src = MeshInfo(X, T, 200);
fprintf('done (found %d vertices)\n',Src.nv);toc;

fprintf('reading the target shape...');tic;
[X, T] = readOff('./Mesh/tr_reg_090');
Tar = MeshInfo(X, T, 200);
fprintf('done (found %d vertices)\n',Tar.nv);toc;

% a few landmark correspondences (to avoid symmetry flipping).
landmarks1 = (500:1000:6000)';
landmarks2 = landmarks1;

landmarks = [landmarks1 landmarks2(:,1)];

SrcLaplaceBasis = Src.laplaceBasis; SrcEigenvalues = Src.eigenvalues;
TarLaplaceBasis = Tar.laplaceBasis; TarEigenvalues = Tar.eigenvalues;
Src.laplaceBasis = SrcLaplaceBasis(:,1:numEigsSrc); Src.eigenvalues = SrcEigenvalues(1:numEigsSrc);
Tar.laplaceBasis = TarLaplaceBasis(:,1:numEigsTar); Tar.eigenvalues = TarEigenvalues(1:numEigsTar);

%% Descriptors
fct_src = [];
fprintf('Computing the descriptors...\n');tic;
fct_src = [fct_src, waveKernelSignature(SrcLaplaceBasis, SrcEigenvalues, Src.Ae, 200)];
fct_src = [fct_src, waveKernelMap(SrcLaplaceBasis, SrcEigenvalues, Src.Ae, 200, landmarks(:,1))];

fct_tar = [];
fct_tar = [fct_tar, waveKernelSignature(TarLaplaceBasis, TarEigenvalues, Tar.Ae, 200)];
fct_tar = [fct_tar, waveKernelMap(TarLaplaceBasis, TarEigenvalues, Tar.Ae, 200, landmarks(:,2))];

% Subsample descriptors (for faster computation). More descriptors is
% usually better, but can be slower. 
fct_src = fct_src(:,1:10:end);
fct_tar = fct_tar(:,1:10:end);

fprintf('done computing descriptors (%d on source and %d on target)\n',size(fct_src,2),size(fct_tar,2)); toc;

assert(size(fct_src,2)==size(fct_tar,2));

% Normalization
no = sqrt(diag(fct_src'*Src.Ae*fct_src))';
fct_src = fct_src ./ repmat(no, [Src.nv,1]);
fct_tar = fct_tar ./ repmat(no, [Tar.nv,1]);

fprintf('Pre-computing the multiplication operators...');tic;
%% Multiplication Operators
numFct = size(fct_src,2);
OpSrc = cell(numFct,1);
OpTar = cell(numFct,1);
for i = 1:numFct
    OpSrc{i} = Src.laplaceBasis'*Src.Ae*(repmat(fct_src(:,i), [1,numEigsSrc]).*Src.laplaceBasis);
    OpTar{i} = Tar.laplaceBasis'*Tar.Ae*(repmat(fct_tar(:,i), [1,numEigsTar]).*Tar.laplaceBasis);
end

Fct_src = Src.laplaceBasis'*Src.Ae*fct_src;
Fct_tar = Tar.laplaceBasis'*Tar.Ae*fct_tar;
fprintf('done\n');toc;

%% Fmap Computation
fprintf('Optimizing the functional map...\n');tic;
Dlb = (repmat(Src.eigenvalues, [1,numEigsTar]) - repmat(Tar.eigenvalues', [numEigsSrc,1])).^2;
Dlb = Dlb/norm(Dlb, 'fro')^2;
constFct = sign(Src.laplaceBasis(1,1)*Tar.laplaceBasis(1,1))*[sqrt(sum(Tar.area)/sum(Src.area)); zeros(numEigsTar-1,1)];

a = 1e-1; % Descriptors preservation
b = 1;    % Commutativity with descriptors
c = 1e-3; % Commutativity with Laplacian 
funObj = @(F) deal( a*sum(sum((reshape(F, [numEigsTar,numEigsSrc])*Fct_src - Fct_tar).^2))/2 + b*sum(cell2mat(cellfun(@(X,Y) sum(sum((X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y).^2)), OpTar', OpSrc', 'UniformOutput', false)), 2)/2 + c*sum( (F.^2 .* Dlb(:))/2 ),...
            a*vec((reshape(F, [numEigsTar,numEigsSrc])*Fct_src - Fct_tar)*Fct_src') + b*sum(cell2mat(cellfun(@(X,Y) vec(X'*(X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y) - (X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y)*Y'), OpTar', OpSrc', 'UniformOutput', false)), 2) + c*F.*Dlb(:));
funProj = @(F) [constFct; F(numEigsTar+1:end)];

F_lb = zeros(numEigsTar*numEigsSrc, 1); F_lb(1) = constFct(1);

% Compute the optional functional map using a quasi-Newton method.
options.verbose = 1;
F_lb = reshape(minConf_PQN(funObj, F_lb, funProj, options), [numEigsTar,numEigsSrc]);
fprintf('done fmap optimization.\n');toc;

%%
fprintf('ICP refinement...');tic;
[F_lb2, ~] = icp_refine(Src.laplaceBasis, Tar.laplaceBasis, F_lb, 5);
fprintf('done\n');toc;

%% Evaluation
% Compute the p2p map

fprintf('Converting to p2p map...');tic;
% fmap before ICP (for comparison)
pF_lb = knnsearch((F_lb*Src.laplaceBasis')', Tar.laplaceBasis);
% fmap after ICP 
pF_lb2 = knnsearch((F_lb2*Src.laplaceBasis')', Tar.laplaceBasis);
fprintf('done\n');toc;

%%
% Plot the results for a random smooth function

fprintf('Visualizing results\n');

figure(1);
nbasis = size(Src.laplaceBasis,2);
kf = rand(nbasis,1)-0.5;
kf = kf.*((nbasis:-1:1).^2');

f = Src.laplaceBasis*kf;
subplot(1,4,1);
trimesh(Src.T, Src.X(:,1), Src.X(:,2), Src.X(:,3),f,'FaceColor','interp');
view([0 90]);
title('f = random smooth function');
axis equal;
axis off;

subplot(1,4,2);
trimesh(Tar.T, Tar.X(:,1), Tar.X(:,2), Tar.X(:,3),f(pF_lb),'FaceColor','interp');
view([0 90]);
title('f(p2pmap) before ICP');
axis equal;
axis off;

subplot(1,4,3);
trimesh(Tar.T, Tar.X(:,1), Tar.X(:,2), Tar.X(:,3),f(pF_lb2),'FaceColor','interp');
view([0 90]);
title('f(p2pmap) after ICP');
axis equal;
axis off;

subplot(1,4,4);
trimesh(Tar.T, Tar.X(:,1), Tar.X(:,2), Tar.X(:,3),Tar.laplaceBasis*F_lb2*kf,'FaceColor','interp');
view([0 90]);
title('Fmap*f (without p2p)');
axis equal;
axis off;

%%
figure(2);
subplot(1,2,1);
samples = euclidean_fps(Src,200);
visualize_map_lines(Tar, Src, pF_lb, samples);
view([0 90]);
title('Map before ICP');
axis equal;
axis off;
fprintf('Mean Euclidean map error (without ICP): %f\n', mean(sqrt(sum((Tar.X-Tar.X(pF_lb,:)).^2,2)))/Tar.sqrt_area);

subplot(1,2,2);
samples = euclidean_fps(Src,200);
visualize_map_lines(Tar, Src, pF_lb2, samples);
view([0 90]);
title('Map after ICP');
axis equal;
axis off;
fprintf('Mean Euclidean map error (with ICP): %f\n', mean(sqrt(sum((Tar.X-Tar.X(pF_lb2,:)).^2,2)))/Tar.sqrt_area);

