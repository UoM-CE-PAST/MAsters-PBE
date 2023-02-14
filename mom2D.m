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
%
% Purpose: Defines a function handle containing the ode system
% corresponding to the two-dimensional method of moments.
% References: 
% (1) Ramkrishna, D., 2000. Population balances : theory and applications
% to particulate systems in engineering. Academic Press.
%
% Input arguments:
% - y: Column vector containing all the dependent variables of the ode
% system
% - solubility: Scalar representing the solubility at a particular
% temperature
% - k11: Parameter used in defining the growth rate perpendicular to the
% radius, G1
% - k12: Parameter used in defining the growth rate perpendicular to the
% radius, G1
% - k21: Parameter used in defining the growth rate perpendicular to the
% length, G2
% - k22: Parameter used in defiing the growth rate perpendicular to the
% length, G2
% - shapeFactor: Factor that relates the m21 cross-moment with the particle
% volume
% - particleDensity: Scalar representing the density of the crystal
% particles
%
% Output arguments:
% - ddt: system of ODEs corresponding to the two-dimensional method of
% moments
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddt = mom2D(y,solubility,k11,k12,k21,k22,shapeFactor,particleDensity)
%% Store variables to solve

m00 = y(1);
m10 = y(2);
m01 = y(3);
m20 = y(4);
m11 = y(5);
m21 = y(6);
c = y(7);

%% Necessary parameters

S = c/solubility;
G1 = k11*(S-1)^k12;
G2 = k21*(S-1)^k22;

%% ODE system

dm00dt = 0; %for the case of no nucleation, this shouldn't change - used to verify simulation
dm10dt = G1*m00;
dm01dt = G2*m11;
dm20dt = 2*G1*m11;
dm11dt = G1*m01 + G2*m10;
dm21dt = 2*G1*m11 + G2*m20;
dcdt = -particleDensity*shapeFactor*dm21dt;

ddt = [dm00dt; dm10dt; dm01dt; dm20dt; dm11dt; dm21dt; dcdt];
end