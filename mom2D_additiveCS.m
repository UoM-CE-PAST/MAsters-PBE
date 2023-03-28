%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
%
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
%
% Last modified:
% - 2023/03/27, MA: Initial creation
%
% Purpose: Defines a function handle containing the ode system
% corresponding to the two-dimensional method of moments in the presence of
% additives. This version is designed to implement a constant
% supersaturation mode of operation.
% References:
% (1) Ramkrishna, D., 2000. Population balances : theory and applications
% to particulate systems in engineering. Academic Press.
% (2) Vetter, T., Mazzotti, M., Brozio, J., 2011. Slowing the growth rate
% of ibuprofen crystals using the polymeric additive pluronic F127. Crystal
% Growth and Design 11. https://doi.org/10.1021/cg200352u
% (3) Salvatori, F., Binel, P., Mazzotti, M., 2019. Efficient assessment of
% combined crystallization, milling, and dissolution cycles for crystal
% size and shape manipulation. Chemical Engineering Science: X 1.
% https://doi.org/10.1016/j.cesx.2018.100004
%
% Input arguments:
%
%
% Output arguments:
% - ddt: system of ODEs corresponding to the two-dimensional method of
% moments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ddt, temperature_mom] = mom2D_additiveCS(t_mom,y,kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,shapeFactor,particleDensity,initialConcentration,initialm21,yieldFactor,maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor)
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

% limits
Tmax = (1/0.036)*log(initialConcentration/(3.37*solubilityFactor));
Smin = supersaturationLimits(1);
Smax = supersaturationLimits(2);
m21min = initialm21;
m21max = yieldFactor*m21min;
cmax = initialConcentration;
cmin = cmax - particleDensity*shapeFactor*(m21max - m21min);
Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

persistent stage
persistent cycleNumber
if isempty(cycleNumber)
    cycleNumber = -1;
end

% temperature calculation dependent on limits:
% switch to a growth stage if minimum volume reached
if m21 <= m21min
    if nargout>1; cycleNumber = cycleNumber + 1; end
    temperature = (1/0.036)*log(concentration_mom/(3.37*solubilityFactor*Smax));
    if temperature < Tmin % enforce lower temperature limit
        temperature = Tmin;
    end
    stage = 1; % growth stage
end

% switch to a dissolution stage if max volume exceeded:
if concentration_mom <= 1.1*cmin && cycleNumber < maxCycleNumber
    temperature = (1/0.036)*log(concentration_mom/(3.37*solubilityFactor*Smin));
    if temperature > Tmax % enforce upper temperature limit
        temperature = Tmax;
    end
    stage = 2; % dissolution stage
end

%normal operation within limits
if stage == 1
    temperature = (1/0.036)*log(concentration_mom/(3.37*solubilityFactor*Smax));
    if temperature < Tmin % enforce lower temperature limit
        temperature = Tmin;
    end
elseif stage == 2
    temperature = (1/0.036)*log(concentration_mom/(3.37*solubilityFactor*Smin));
    if temperature > Tmax % enforce upper temperature limit
        temperature = Tmax;
    end
end

% solubility and supersaturation calculation
solubility = solubilityFactor*3.37*exp(0.036*temperature);
supersaturation = concentration_mom/solubility;
% growth/dissolution rate dependent on supersaturation
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

dm00dt = 0; %for the case of no nucleation, this shouldn't change
% - used to verify simulation

% other moments:
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

dcdt = -particleDensity*shapeFactor*dm21dt; % concentration

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

% used to store the temperature used by ode15s
if nargout>1; temperature_mom = temperature; end
end