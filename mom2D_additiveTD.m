%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
% 
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
% 
% Last modified:
% - 2023/03/20, MA: Initial creation
%
% Purpose: Defines a function handle containing the ode system
% corresponding to the two-dimensional method of moments in the presence of
% additives.
% References: 
% (1) Ramkrishna, D., 2000. Population balances : theory and applications
% to particulate systems in engineering. Academic Press.
% (2) Vetter, T., Mazzotti, M., Brozio, J., 2011. Slowing the growth rate
% of ibuprofen crystals using the polymeric additive pluronic F127. Crystal
% Growth and Design 11. https://doi.org/10.1021/cg200352u
%
% Input arguments:
%
%
% Output arguments:
% - ddt: system of ODEs corresponding to the two-dimensional method of
% moments
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddt = mom2D_additiveTD(t_mom,y,kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,temperatureRamp,shapeFactor,particleDensity,solubilityFactor,growthFactor)
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

m03 = y(12);
m40 = y(13);
m13 = y(14);
m23 = y(15);
m41 = y(16);

concentration_mom = y(17);

%% Necessary parameters
temperature = interp1(temperatureRamp(1,:),temperatureRamp(2,:),t_mom);  
solubility = solubilityFactor*3.37*exp(0.036*temperature);
supersaturation = concentration_mom/solubility;
if supersaturation > 1
    G1 = kg11*exp(-kg12/(temperature+273.15))*(supersaturation-1)^kg13;
    G2 = growthFactor*kg21*exp(-kg22/(temperature+273.15))*(supersaturation-1)^kg23;
elseif supersaturation < 1
    G1 = kd11*exp(-kd12/(temperature+273.15))*(supersaturation-1);
    G2 = kd21*exp(-kd22/(temperature+273.15))*(supersaturation-1);
else
    G1 = 0;
    G2 = 0;
end

%% ODE system

dm00dt = 0; %for the case of no nucleation, this shouldn't change - used to verify simulation
dm10dt = G1*m00;
dm01dt = G2*m00;
dm20dt = 2*G1*m10;
dm02dt = 2*G2*m01;
dm11dt = G1*m01 + G2*m10;
dm21dt = 2*G1*m11 + G2*m20;
dm12dt = G1*m02 + 2*G2*m11;
dm22dt = 2*G1*m12 + 2*G2*m21;
dm30dt = 3*G1*m20;
dm31dt = 3*G1*m21 + G2*m30;
dm03dt = 3*G2*m02;
dm40dt = 4*G1*m30;
dm13dt = G1*m03 + 3*G2*m12;
dm23dt = 2*G1*m13 + 3*G2*m22;
dm41dt = 4*G1*m31 + G2*m40;

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
    dm03dt;
    dm40dt;
    dm13dt;
    dm23dt;
    dm41dt;
    dcdt];
end