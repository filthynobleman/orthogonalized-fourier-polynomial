% Plots the first eigenfunctions, eigenproducts and orthogonalized
% eigenproducts for the real line.
%
% This script reproduces Fig. 7 of the reference paper.
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
MaxK = 4;
syms x;
PlotAnalytic = false;
NumPlotPoints = 50;


%% Create the functions
% Eigenfunctions
K = (0:MaxK);
Mu = K;
Phi = sin(K .* x);
% Eigenproducts
PolyPhi = [x];
for i = 2:length(K)
    for j = i:length(K)
        PolyPhi(end + 1) = Phi(i) * Phi(j);
    end
end
PolyPhi = PolyPhi(2:end);
% Orthogonal
Ortho = PolyPhi;
for i = 2:length(Ortho)
    f = Ortho(i);
    Coeff = zeros(length(K) + i - 1, 1);
    for j = 1:length(K)
        Coeff(j) = double(vpaintegral(f * Phi(j), [0, 2 * pi]));
    end
    for j = 1:i-1
        Coeff(length(K) + j) = double(vpaintegral(f * Ortho(j), [0, 2 * pi]));
    end
    Ortho(i) = f - [Phi, Ortho(1:i-1)] * Coeff;
end
for i = 1:length(Ortho)
    Norm = double(vpaintegral(Ortho(i)^2, [0, 2 * pi]));
    Ortho(i) = Ortho(i) / sqrt(Norm);
end






%% Plot analytic
if PlotAnalytic
    % Plot eigenfunctions
    for i = 1:length(Phi)
        subplot(5, length(K), i);
        fplot(Phi(i), [0, 2 * pi], 'b', ...
              'LineWidth', 2);
        ylim([-1, 1]);
    end

    % Plot eigenproducts
    for i = 1:length(PolyPhi)
        subplot(5, length(K), i + length(K));
        fplot(PolyPhi(i), [0, 2 * pi], 'r', ...
              'LineWidth', 2);
        ylim([-1, 1]);
    end

    % Plot eigenproducts
    for i = 1:length(Ortho)
        subplot(5, length(K), i + 3 * length(K));
        fplot(Ortho(i), [0, 2 * pi], 'g', ...
              'LineWidth', 2);
        ylim([-1, 1]);
    end
else
    X = linspace(0, 2 * pi, NumPlotPoints);
    
    % Plot eigenfunctions
    for i = 1:length(Phi)
        subplot(5, length(K), i);
        plot(X, double(subs(Phi(i), x, X)), 'b', ...
              'LineWidth', 2);
        ylim([-1, 1]);
        xlim([0, 2*pi]);
    end

    % Plot eigenproducts
    for i = 1:length(PolyPhi)
        subplot(5, length(K), i + length(K));
        plot(X, double(subs(PolyPhi(i), x, X)), 'r', ...
              'LineWidth', 2);
        ylim([-1, 1]);
        xlim([0, 2*pi]);
    end

    % Plot eigenproducts
    for i = 1:length(Ortho)
        subplot(5, length(K), i + 3 * length(K));
        plot(X, double(subs(Ortho(i), x, X)), 'g', ...
              'LineWidth', 2);
        ylim([-1, 1]);
        xlim([0, 2*pi]);
    end
end

















































