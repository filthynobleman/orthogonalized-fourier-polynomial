function [fs,idxs_func] = funcset(M, k, K, diam)



fs = [];
idxs_func =[0];

% hk k (1)
x = randi(M.n);
t = 1e-3;
hk = M.Phi(:, 1:k) * diag(exp(- M.Lambda(1:k) * t)) * M.Phi(:, 1:k)';
add = hk(x, :)';
add = add ./ max(abs(add));
fs = [fs, add]; clear add;
idxs_func = [idxs_func, size(fs,2)];

% hk K (2)
x = randi(M.n);
t = 1e-3;
hk = M.Phi(:, 1:K) * diag(exp(- M.Lambda(1:K) * t)) * M.Phi(:, 1:K)';
add = hk(x,:)';
add = add ./ max(abs(add));
fs = [fs, add]; clear add;
idxs_func = [idxs_func, size(fs, 2)];

% HKS (3) 50
tmin = abs(4*log(10) / M.Lambda(k + 1));
tmax = abs(4*log(10) / M.Lambda(2));
nstep = 99;
stepsize = (log(tmax) - log(tmin)) / nstep;
logts = log(tmin):stepsize:log(tmax);
t = exp(logts);
HKS = zeros(M.n, length(t));
% skipping the first freq. as the authors do in their code
for i=1:length(t)
    HKS(:,i) = sum(...
        (M.Phi(:,2:k+1)).^2 .* repmat(exp(-t(i) * M.Lambda(2:k+1))', M.n, 1), ...
        2);
end
HKS_idxs = randperm(50);
add =  HKS(:, HKS_idxs(1:10));
add = add ./ max(abs(add));
fs = [fs, add]; clear add;
idxs_func = [idxs_func, size(fs, 2)];

% WKS (4) 20
WKS = waveKernelSignature_exact(M.Phi(:, 1:k+1), M.Lambda(1:k+1), 100);
WKS_idxs = randperm(100);
add =  WKS(:, WKS_idxs(1:10));
add = add ./ max(abs(add));
fs = [fs, add]; clear add;
idxs_func = [idxs_func, size(fs,2)];

% Random (5)
add = rand(M.n, 1);
add = add ./ max(abs(add));
fs = [fs, add]; clear add;
idxs_func = [idxs_func, size(fs,2)];

% coordinates (6) 3
add =  M.X;
add = add ./ max(abs(add));
fs = [fs, add]; clear add;
idxs_func = [idxs_func, size(fs,2)];

% indicator (7) 1
nset = 10;
ind = zeros(nset, 1);
vor = ones(M.n, 1);
% geo = dijkstra_to_all(M, [1]);
geo = mesh.metrics.dijkstra(M, 1);
[~, index] = max(geo);
ind(1) = index;
% geo = dijkstra_to_all(M, [index]);
geo = mesh.metrics.dijkstra(M, index);
for vv = 1:(nset-1)
    [~, index] = max(geo);
    ind(vv + 1) = index;
%     local_geo = dijkstra_to_all(M, ind(vv + 1));
    local_geo = mesh.metrics.dijkstra(M, ind(vv + 1));
    vor(local_geo < geo) = (vv + 1);
    geo = min(geo, local_geo);
end
i_vor = randi([1, nset]); 
f = (vor == i_vor);
add =  f;
add = add ./ max(abs(add));
fs = [fs, add]; clear add;
idxs_func = [idxs_func, size(fs, 2)];

% shot (8)
r = 0.04*diam;
d = calc_shot([M.X M.Y M.Z]',M.TRIV', [1:length(M.X)], 9, r, 3);
shot_idxs = randperm(size(d' ,2));
add =  d(shot_idxs(1:10),:)';
add = add./max(abs(add));
fs = [fs, add]; clear add;
idxs_func = [idxs_func, size(fs,2)];

% AWFT
AWFTparam.freqs0 = [1 1];
AWFTparam.taus = [0.002 0.05]; 
AWFTparam.angles = [-1 90 180]; % -1 is for the isotropic case
AWFTparam.curv_smooth = 10;
AWFTparam.n_eigen = 100;
AWFTparam.alpha = [ 100];
AWFTfunctions{1} = 'ShapeIndex';
AWFTfunctions{2} = 'Fiedler';
AWFT = compute_AWFTdesc(M,AWFTfunctions,AWFTparam);
add =  AWFT;
add = add./max(abs(add));
fs = [fs, add]; clear add;
idxs_func = [idxs_func, size(fs,2)];

% % hk K (2)
% for jjj = 1:10 
%     x = randi(M.n);
%     t = 1e-6;
%     hk = M.Phi * diag(exp(- M.Lambda * t)) * M.Phi';
%     add = hk(x, :)';
%     add = add ./ max(abs(add));
%     fs = [fs, add]; clear add;
% end
% idxs_func = [idxs_func, size(fs, 2)];
end


function wks = waveKernelSignature_exact(laplaceBasis, eigenvalues, numTimes)
% This method computes the wave kernel signature for each vertex on a list.
% It uses precomputed LB eigenstuff stored in "mesh" and automatically
% chooses the time steps based on mesh geometry.

numEigenfunctions = size(eigenvalues,1);

% D = laplaceBasis' * (Ae * laplaceBasis.^2);
D = laplaceBasis.^2;

absoluteEigenvalues = abs(eigenvalues);
emin = log(absoluteEigenvalues(2));
emax = log(absoluteEigenvalues(end));
s = 7*(emax-emin) / numTimes; % Why 7?
emin = emin + 2*s;
emax = emax - 2*s;
es = linspace(emin,emax,numTimes);

T = exp(-(repmat(log(absoluteEigenvalues),1,numTimes) - ...
    repmat(es,numEigenfunctions,1)).^2/(2*s^2));
wks = D*T;
% wks = laplaceBasis*wks;
end