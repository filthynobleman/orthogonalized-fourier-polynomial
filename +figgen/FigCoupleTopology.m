% Plots two figures, highlighting their topology.
%
% This script reproduces Fig. 11 of the reference paper.
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
MeshFiles = ["Meshes/visible_mesh/tr_reg_065_remeshed5K.off";
             "Meshes/visible_mesh/tr_reg_081.off"];
         

%% Load the meshes
fprintf("Loading the meshes... ");
tic;
Src = mesh.init(MeshFiles(1));
Src = mesh.transform.rotate(Src, 'x', 90);
Trg = mesh.init(MeshFiles(2));
Trg = mesh.transform.rotate(Trg, 'x', 90);
toc;



%% Plot
subplot(1, 2, 1);
utils.plot_scalar_map(Src, zeros(Src.n, 1), white, true, 'faceted', false);
subplot(1, 2, 2);
utils.plot_scalar_map(Trg, zeros(Trg.n, 1), white, true, 'faceted', false);













































































