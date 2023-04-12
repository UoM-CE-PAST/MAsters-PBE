%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
% 
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
% 
% Last modified:
% - 2022/03/06, MA: initial creation
% - 2023/03/07, MA: Fixed growth rate
% - 2023/04/30, MA: changed kinetics + added average length support
%
% Purpose: Defines a function handle containing the ode system
% corresponding to the one-dimensional method of moments.This function
% differentiates itself by using the solubility and growth rate outlined in
% the paper by T. Vetter et al.
% References: 
% (1) Ramkrishna, D., 2000. Population balances : theory and applications
% to particulate systems in engineering. Academic Press.
% (2) Vetter, T., Mazzotti, M., Brozio, J., 2011. Slowing the growth rate
% of ibuprofen crystals using the polymeric additive pluronic F127. Crystal
% Growth and Design 11. https://doi.org/10.1021/cg200352u
% (3) Vetter, T., Mazzotti, M., Brozio, J., 2011. Slowing the growth rate
% of ibuprofen crystals using the polymeric additive pluronic F127. Crystal
% Growth and Design 11. https://doi.org/10.1021/cg200352u
%
% Input arguments:
% - y: Column vector containing all the dependent variables of the ode
% system
% - solubility: Scalar representing the solubility at a particular
% temperature
% - k1: Scalar reperesenting one of the growth rate parameters [m/s] (3).
% - k2: Scalar representing another one of the growth rate parameters [K]
% (3).
% - k3: Scalar representing another one of the growth rate parameter [K^2]
% (3).
% p0, p1, p2, p3, p4: scalars representing the parameters needed for the
% solubility polynomial [-] (4).
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
function ddt = mom_additive(t_mom, y, temperatureRamp, kg1, kg2, kg3, kd1, kd2, particleDensity, ...
    shapeFactor, growthFactor, solubilityFactor)
%% Store variables to be solved

m0 = y(1);
m1 = y(2);
m2 = y(3);
m3 = y(4);
m4 = y(5);
m5 = y(6);
concentration_mom = y(7);

%% Necessary parameters
temperature = interp1(temperatureRamp(1,:),temperatureRamp(2,:),t_mom);
solubility = solubilityFactor*3.37*exp(0.036*temperature);
supersaturation = concentration_mom/solubility;
if supersaturation > 1
    G = growthFactor*kg1*exp(-kg2/(temperature+273.15))*(supersaturation-1)^kg3; %[um/h]
elseif supersaturation < 1
    G = kd1*exp(-kd2/(temperature+273.15))*(supersaturation-1); %[um/h]
else
    G = 0;
end

%% ODE system

dm0dt = 0; % This is added to double check the simulation is correct
dm1dt = G*m0;
dm2dt = 2*G*m1;
dm3dt = 3*G*m2;
dm4dt = 4*G*m3;
dm5dt = 5*G*m4;
dcdt = -particleDensity*shapeFactor*dm3dt;

ddt = [dm0dt; dm1dt; dm2dt; dm3dt; dm4dt; dm5dt; dcdt];
end