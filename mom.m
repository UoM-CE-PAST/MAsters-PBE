%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
% 
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
% 
% Last modified:
% - 2022/10/10, MA: initial creation
% - 2022/02/23, MA: improved readability
% - 2022/02/24, MA: added temperature dependence
% - 2022/02/24, MA: added constant supersaturation mode
%
% Purpose: Defines a function handle containing the ode system
% corresponding to the one-dimensional method of moments.
% References: 
% (1) Ramkrishna, D., 2000. Population balances : theory and applications
% to particulate systems in engineering. Academic Press.
%
% Input arguments:
% - y: Column vector containing all the dependent variables of the ode
% system
% - solubility: Scalar representing the solubility at a particular
% temperature
% - k1: growth rate parameter
% - k2: growth rate parameter
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
function ddt = mom(t_mom, y, temperatureRamp, k1, k2, particleDensity, ...
    shapeFactor, operationMode)
%% Store variables to be solved

m0 = y(1);
m1 = y(2);
m2 = y(3);
m3 = y(4);
concentration_mom = y(5);

%% Necessary parameters

if operationMode == 1 % time dependent operation
    temperature = interp1(temperatureRamp(1,:),temperatureRamp(2,:),t_mom);  
elseif operationMode == 2 % concentration dependent operation
    temperature = interp1(temperatureRamp(1,:),temperatureRamp(2,:),concentration_mom);
end
solubility = 3.37*exp(0.0359*temperature);
S = concentration_mom/solubility;
G = k1*(S-1)^k2;
%% ODE system

dm0dt = 0; % This is added to double check the simulation is correct
dm1dt = G*m0;
dm2dt = 2*G*m1;
dm3dt = 3*G*m2;
dcdt = -particleDensity*shapeFactor*dm3dt;

ddt = [dm0dt; dm1dt; dm2dt; dm3dt; dcdt];
end