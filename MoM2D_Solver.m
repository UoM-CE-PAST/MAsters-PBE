%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
% 
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
% 
% Last modified:
% - 2023/02/13, MA: Initial creation
% - 2023/02/14, MA: minor fixes
%
% Purpose: Defines a function handle containing the ode system
% corresponding to the two-dimensional method of moments.
% References: 
% (1) Ramkrishna, D., 2000. Population balances : theory and applications
% to particulate systems in engineering. Academic Press.
%
% Input arguments
% 
% Output arguments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc
close all

%% Define initial PSD
dL1 = 1;
L1 = 1:dL1:1000;
dL2 = 1;
L2 = 1:dL2:1000;
[l1,l2] = meshgrid(L1,L2);

L = [l1(:) l2(:)];
mean = [200 400];
standardDeviation = [20 0; 0 20];

f0 = 1e5*mvnpdf(L,mean, standardDeviation);
f0 = reshape(f0,length(L2),length(L1));

surf(L1,L2,f0)
clim([min(f0(:))-0.5*range(f0(:)),max(f0(:))])
axis([180 220 380 420 0 800])

%% Calculate initial cross-moments (using integration)

initialm00 = trapz(L2,trapz(L1,f0));
initialm10 = trapz(L2,trapz(L1,L2'.*f0));
initialm01 = trapz(L2,trapz(L1,L1.*f0));
initialm20 = trapz(L2,trapz(L1.^2.*f0));
initialm11 = trapz(L2,trapz(L1,L1.*L2'.*f0));
initialm21 = trapz(L2, trapz(L1,L1.^2.*L2'.*f0));

%% Constant values

T = 25;
solubility = 3.37*exp(0.0359*T);
c0 = 30;
k11 = 5;
k12 = 1;
k21 = 2.5;
k22 = 1;
shapeFactor = pi();
particleDensity = 1.46e-12;

%% ODE system integration

%Initial conditions
y0 = [initialm00 initialm10 initialm01 initialm20 initialm11 initialm21 c0];

timeSpan = [0 10];

[t, y] = ode45(@(t,y)mom2D(y,solubility,k11,k12,k21,k22,shapeFactor,particleDensity), timeSpan, y0);

figure(2)
plot(t,y(:,7))