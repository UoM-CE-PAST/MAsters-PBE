%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
% 
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
% 
% Last modified:
% - 2023/02/10, MA: Initial creation
% - 2023/02/13, MA: Added ability to model dissolution
%
% Purpose: This file is used to track/visualise changes in the particle
% size distribution over a given time range as simulated by the function:
% HighRes1D.
%
% References: 
% (1) LeVeque, R.J., 2002. Finite Volume Methods for Hyperbolic Problems, Cambridge Texts in Applied Mathematics. Cambridge University Press, Cambridge. https://doi.org/10.1017/CBO9780511791253
% (2) Gunawan, R., Fusman, I., Braatz, R.D., 2004. High resolution algorithms for multidimensional population balance equations. AIChE Journal 50. https://doi.org/10.1002/aic.10228
% (3) Ma, D.L., Tafti, D.K., Braatz, R.D., 2002. High-resolution simulation of multidimensional crystal growth. Industrial and Engineering Chemistry Research 41. https://doi.org/10.1021/ie010680u
%
% Input Arguments:
% dL: Scalar representing the length of the length step
%
% L: 1d array representing the spatial domain 
%
% tmax: Scalar representing the duration of the simulation
% 
% k1: Scalar reperesenting one of the growth rate parameters
%
% k2: Scalar representing another one of the growth rate parameters
%
% shapeFactor: Scalar representing particle shape factor
%
% temperature: Scalar representing the temperature
%
% ParticleDensity: Scalar representing the particle density
%
% c0: Scalar representing the initial concentration
%
% f0: Scalar representing the initial particle distribution
%
% Output arguments:
% f: 2d array containing the particle size distribution at every time and
% length.
%
% c: 1d array containing the concentration of the liquid phase at each
% time step
%
% G: 1d array containing the particle growth rate at each time step
%
% S: 1d array containing the supersaturation of the liquid phase at each
% time step
%
% m3: 1d array copntaing the 3rd moment of the particle distribution
% (proportional to particle volume)
%
% t: 1d array containing the time elapsed since the start of the
% simulation for each time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc

%% Essential parameters
temperature=25; % C
c0=50; % g/kg
k1 = 10; % um/h
k2 = 1;
shapeFactor = 1; % assuming cuboidal particles
particleDensity = 1.46e-12; % g/um3
tmax=100; %h range of t required

%Length range and length step
dL=1; %um
L=1:dL:1000; %um

%Initial PSD

% Gaussian
f0 = 1e5*normpdf(L,300,20); 

% Pulse
% f0 = zeros(1,length(L));
% f0(100:200)=1e3; 

%% Base Case

[f, c, G, S, m3, t, ce] = HighRes1D(dL, L, tmax, k1, k2, shapeFactor, temperature, particleDensity, c0, f0);

%Static plot

plot(L,f(:,1), 'linewidth',1.2), xlabel('Length [μm]'), ylabel('f [μm^{-1} kg^{-1}]'), hold on, plot(L,f(:,end), 'linewidth',1.2), xlabel('Length [μm]'), ylabel('f [μm^{-1} kg^{-1}]')
legend('Inital PSD','Final PSD')
set(gca,'FontSize',18)
title('(a)','FontSize',24)