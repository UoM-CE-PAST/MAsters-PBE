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
% - 2023/03/13, MA: added temperature dependence
% - 2023/03/14, MA: major fixes to ODE system
% - 2023/03/21, MA: added support for volume-weighted average
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
function ddt = mom2D(t_mom,y,k11,k12,k21,k22,temperatureRamp,shapeFactor,particleDensity,operationMode)
%% Store variables to solve

m00 = y(1);
m10 = y(2);
m01 = y(3);
m20 = y(4);
m02 = y(5);
m11 = y(6);
m21 = y(7);
m12 = y(8);
m22 = y(9);
m30 = y(10);
m31 = y(11);
concentration_mom = y(12);

%% Necessary parameters

if operationMode == 1 % time dependent operation
    temperature = interp1(temperatureRamp(1,:),temperatureRamp(2,:),t_mom);  
elseif operationMode == 2 % concentration dependent operation
    temperature = interp1(temperatureRamp(1,:),temperatureRamp(2,:),concentration_mom);
end
solubility = 3.37*exp(0.0359*temperature);
supersaturation = concentration_mom/solubility;
G1 = k11*(supersaturation-1)^k12;
G2 = k21*(supersaturation-1)^k22;

%% ODE system

dm00dt = 0; %for the case of no nucleation, this shouldn't change - used to verify simulation
dm10dt = G1*m00;
dm01dt = G2*m00;
dm20dt = 2*G1*m10;
dm02dt = 2*G2*m01;
dm11dt = G1*m01 + G2*m10;
dm21dt = 2*G1*m11 + G2*m20;
dm12dt = G1*m02 + 2*G2m11;
dm22dt = 2*G1*m12 + 2*G2*m21;
dm30dt = 3*G1*m20;
dm31dt = 3*G1*m21 + G2*m30;

dcdt = -particleDensity*shapeFactor*dm21dt;

ddt = [dm00dt;
    dm10dt;
    dm01dt;
    dm20dt;
    dm02dt;
    dm11dt;
    dm21dt;
    dm12dt;
    dm22dt;
    dm30dt;
    dm31dt;
    dcdt];
end