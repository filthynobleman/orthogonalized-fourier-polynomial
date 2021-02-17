% Shape Analysis with Anisotropic Windowed Fourier Transform (3DV 2016)
% CODE for the computation of WFT and AWFT descriptors
% AUTHORS:  Simone  Melzi , simone.melzi@univr.it
%           Emanuele  Rodola , emanuele.rodola@usi.ch
%           Umberto  Castellani , umberto.castellani@univr.it
%           Michael  Bronstein , michael.bronstein@usi.ch
% 
% Verona 15 Settembre 2016 
% -------------------------------------------------------------------------
% The output of this code are the WFT and AWFT descriptors presented in the
% paper Shape Analysis with Anisotropic Windowed Fourier Transform
% -------------------------------------------------------------------------
% select input mesh as a struct shape, with the following fields:
% shape.X, shape.Y, shape.Z:  vectors of dimension Nx1 that are the X,Y,Z
%                             coordinates of vertices in R^3
% shape.TRIV: a num_triangles x 3 matrix each row of which is a triangle in
%             the mesh.

clear;
close all;
clc;

% add path of this folder and all its subfolders
addpath(genpath('..\clean_code_AWFT'))

name_shape = 'tr_reg_090'; % 'horse';

% ---------------------------------------------------- editable parameters
numEig = 200; % number of eigenvetors to consider

% WFT parameters
WFTparam.freqs0 = [1 1 1 1 1];
WFTparam.taus = [0.001, 0.005, 0.01, 0.05, 0.1];
WFTparam.n_eigen = numEig;
WFTparam.curv_smooth = 10;

% select set of function to analyze in WFT
WFTfunctions{1} = 'geovec1';
WFTfunctions{2} = 'geovec3';
WFTfunctions{3} = 'geovec4';
WFTfunctions{4} = 'geovec5';
WFTfunctions{5} = 'geovec6';
WFTfunctions{6} = 'Constant';
WFTfunctions{7} = 'Fiedler';
WFTfunctions{8} = 'ShapeIndex';
WFTfunctions{9} = 'HKS';
WFTfunctions{10} = 'WKS';

% AWFT parameters
AWFTparam.freqs0 = [1 1];
AWFTparam.taus = [0.002 0.05]; 
AWFTparam.angles = [-1 45 90 135 180]; % -1 is for the isotropic case
AWFTparam.curv_smooth = 10;
AWFTparam.n_eigen = numEig;
AWFTparam.alpha = [ 100 300 ];

% select set of function to analyze in AWFT
AWFTfunctions{1} = 'geovec3';
AWFTfunctions{2} = 'geovec4';
AWFTfunctions{3} = 'ShapeIndex';
AWFTfunctions{4} = 'Fiedler';
AWFTfunctions{5} = 'Constant';

% ------------------------------------------------------------------------

% dispay infos
fprintf('processing shape %s: \n',name_shape);
fprintf('[i] compute AWFT and WFT descriptors:');

time_start = tic;
    
% load the current shape
load(['.\data\',name_shape,'.mat']);
    
% compute WFT
WFT = compute_WFTdesc(shape,WFTfunctions,WFTparam);
% compute AWFT
AWFT = compute_AWFTdesc(shape,AWFTfunctions,AWFTparam);

% elasped time
elapsed_time = toc(time_start);
fprintf('%3.2fs\n',elapsed_time);