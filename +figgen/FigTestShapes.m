% Plots the four shapes used in the reference paper for testing
% assessements and checking theoretical results.
%
% This script reproduces Fig. 3 of the reference paper.
%
% This code is part of the implementation of the paper "Orthogonalized 
% Fourier Polynomials for Signal Approximation and Transfer", authored by
% Filippo Maggioli, Simone Melzi, Michael Bronstein, Maks Ovsjanikov and
% Emanuele Rodolà.

%% Clear workspace
clear all
close all
clc


%% Parameters
MeshFiles = ["Meshes/tr_reg/tr_reg_000.ply";
             "Meshes/bunny.off";
             "Meshes/cat0.off";
             "Meshes/texturedMeshes/t1_donut.obj"];
         
         
%% Load the shapes
fprintf("Loading the meshes... ");
tic;
M = cell(length(MeshFiles), 1);
for i = 1:length(MeshFiles)
    M{i} = mesh.init(MeshFiles(i));
end
toc;

%% Adjust the shapes orientation
M{3} = mesh.transform.rotate(M{3}, 'x', -90);
M{4} = mesh.transform.rotate(M{4}, 'x', 60);



%% Plot
fig = figure;
fig.WindowState = 'maximized';
for i = 1:4
    subplot(2, 2, i);
    utils.plot_scalar_map(M{i}, zeros(M{i}.n, 1), utils.cmaps.constant(0.9), true);
end























































