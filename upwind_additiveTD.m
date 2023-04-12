%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
%
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
%
% Last modified:
% - 2023/03/06, MA: Initial creation
%
% Purpose: Implements an upwind finite volume method to solve for the
% time evolution of a 1D particle size distributionover a given time range.
% This function differentiates itself by using thesolubility and growth 
% rate outlined in the paper by T. Vetter et al.
%
% References:
% (1) LeVeque, R.J., 2002. Finite Volume Methods for Hyperbolic Problems,
% Cambridge Texts in Applied Mathematics. Cambridge University Press,
% Cambridge. https://doi.org/10.1017/CBO9780511791253
% (2) Gunawan, R., Fusman, I., Braatz, R.D., 2004. High resolution
% algorithms for multidimensional population balance equations. AIChE
% Journal 50. https://doi.org/10.1002/aic.10228
% (3) Ma, D.L., Tafti, D.K., Braatz, R.D., 2002. High-resolution simulation
% of multidimensional crystal growth. Industrial and Engineering Chemistry
% Research 41. https://doi.org/10.1021/ie010680u
% (4) Vetter, T., Mazzotti, M., Brozio, J., 2011. Slowing the growth rate
% of ibuprofen crystals using the polymeric additive pluronic F127. Crystal
% Growth and Design 11. https://doi.org/10.1021/cg200352u
%
% Input Arguments:
% dL: Scalar representing the length of the length step
%
% L: 1d array representing the spatial domain
%
% tmax: Scalar representing the duration of the simulation
%
% k1: Scalar reperesenting one of the growth rate parameters [m/s] (4).
%
% k2: Scalar representing another one of the growth rate parameters [K]
% (4).
%
% k3: Scalar representing another one of the growth rate parameter [K^2]
% (4).
%
% p0, p1, p2, p3, p4: scalars representing the parameters needed for the
% solubility polynomial (4).
%
% kv: Scalar representing particle shape factor
%
% T: Scalar representing the temperature
%
% ParticleDensity: Scalar representing the particle density
%
% c0: Scalar representing the initial concentration
%
% initialPSD: Scalar representing the initial particle distribution
%
% T0: Scalar containing equilibrium temperature for initial concentration
%
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
% m3: 1d array containing the 3rd moment of the particle distribution
% (proportional to particle volume)
%
% t: 1d array containing the time elapsed since the start of the
% simulation for each time step
%
% temperature: 1d array containing the temperature at each time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, concentration, G, supersaturation, m3, t,...
    temperature] = upwind_additiveTD(dL, L, simulationTime, kg1, kg2, kg3, kd1, kd2, shapeFactor,...
    temperatureRamp, particleDensity, initialConcentration, initialPSD, growthFactor, solubilityFactor)

%Initialise time dependent variables
t = zeros(1,5000);
m3 = zeros(1,length(t));
concentration = zeros(1,length(t));
temperature = zeros(1,length(t));
supersaturation = zeros(1,length(t));
G = zeros(1,length(t));
f = zeros(length(L),length(t));


%Initial values
f(:,1)=initialPSD;
m3(1)=trapz(L,L.^3.*f(:,1)');
concentration(1)=initialConcentration;
T0 = (1/0.036)*log(initialConcentration/(3.37*solubilityFactor));
temperature(1) = temperatureRamp(2,1);
solubility = solubilityFactor*3.37*exp(0.036*temperature(1));
supersaturation(1)=initialConcentration/solubility;
if supersaturation(1) > 1
    G(1)=growthFactor*kg1*exp(-kg2/(temperature(1)+273.15))*(supersaturation(1)-1)^kg3; %[um/h]
elseif supersaturation(1) < 1
    G(1) = kd1*exp(-kd2/(temperature(1)+273.15))*(supersaturation(1)-1); %[um/h]
else
    G(1) = 0;
end

%initialise smoothness and flux limiter functions
fluxLimiter = zeros(length(L),1);

%Courant number to specify the maximum stable time step
CourantNumber = 1;

% set reference time and time index
n=1;

% if system starts at equilibrium (i.e. G=0) then ensure system doesn't
% change until driving force is non-zero
while temperature(n) == T0 || G(n) == 0
    t(n+1)=temperatureRamp(1,n+1);
    f(:,n+1)=initialPSD;
    m3(n+1)=trapz(L,L.^3.*f(:,n+1)');
    concentration(n+1)=initialConcentration;
    temperature(n+1) = temperatureRamp(2,n+1);
    solubility = solubilityFactor*3.37*exp(0.036*temperature(n+1));
    supersaturation(n+1)=concentration(n+1)/solubility;
    if supersaturation(n+1) > 1
    G(n+1)=growthFactor*kg1*exp(-kg2/(temperature(n+1)+273.15))*(supersaturation(n+1)-1)^kg3; %[um/h]
    elseif supersaturation(n+1) < 1
        G(n+1) = kd1*exp(-kd2/(temperature(n+1)+273.15))*(supersaturation(n+1)-1); %[um/h]
    else
        G(n+1) = 0;
    end

    n=n+1;
end

%while loop to update and store cell averages over time
if G(n)>0
    % Growth
    while t(n)<simulationTime

        %Calculate stable time step using available values
        dt=CourantNumber*dL/G(n);

        %Check if the max stable time step will exceed time range
        if simulationTime-t(n)<=dt
            t(n+1)=simulationTime;
            dt = simulationTime-t(n);
        else
            t(n+1)=t(n)+dt;
        end

        % Update the PSD at the new time using high resolution method
        %% 1-Inflow boundary

        %Determine smoothness and calculate appropriate flux limiter for cell
        %inlet and outlet fluxes

        f(1,n+1)=f(1,n)-dt/dL*G(n)*f(1,n)-0.5*dt/dL*G(n)*(1-dt/dL*G(n))*(fluxLimiter(2)*(f(2,n)-f(1,n))-fluxLimiter(1)*(f(1,n)));

        %% 2-Interior volume

        f(2:end-1,n+1) = f(2:end-1,n) - dt/dL*(G(n)*f(2:end-1,n)-G(n)*f(1:end-2,n)) - 0.5*dt/dL*G(n)*(1-dt/dL*G(n))*(fluxLimiter(3:end).*(f(3:end,n)-f(2:end-1,n)) - fluxLimiter(2:end-1).*(f(2:end-1,n)-f(1:end-2,n)));

        %% 3-Outflow boundary

        % An outlet flux limiter is not required for the outflow boundary

        f(end,n+1)=f(end,n)-dt/dL*(G(n)*f(end,n)-G(n)*f(end-1,n))-0.5*dt/dL*G(n)*(1-dt/dL*G(n))*(-fluxLimiter(end)*(f(end,n)-f(end-1,n)));

        %% Use liquid phase mass balance to determine supersaturation at next time step
        m3(n+1)=trapz(L,L.^3.*f(:,n+1)');
        concentration(n+1) = concentration(n) - particleDensity*shapeFactor*(m3(n+1)-m3(n));

        % Interpolate temperature to find solubility and superstaturation
        % if temperature is time-dependent, interpolate with time
        temperature(n+1) = interp1(temperatureRamp(1,:),temperatureRamp(2,:),t(n+1));
        solubility = solubilityFactor*3.37*exp(0.036*temperature(n+1));
        supersaturation(n+1)=concentration(n+1)/solubility;

        if supersaturation(n+1)<=1 % Necessary to make sure it remains a growth problem
            supersaturation(n+1)=1;
        end

        G(n+1)=growthFactor*kg1*exp(-kg2/(temperature(n+1)+273.15))*(supersaturation(n+1)-1)^kg3; %[um/h]

        % Increase time counter
        n=n+1;
        
    end

elseif G(n)<0
    % Dissolution
    % Courant number is now negative
    CourantNumber = -CourantNumber;
    while t(n)<simulationTime

        %Calculate stable time step using available values
        dt=CourantNumber*dL/G(n);
   
        %Check if the max stable time step will exceed time range
        if simulationTime-t(n)<=dt
            t(n+1)=simulationTime;
            dt = simulationTime-t(n);
        else
            t(n+1)=t(n)+dt;
        end

        %Calculate the PSD at the new time using high resolution method:
        % for 1, 2 & 3, the smoothness is determined and the appropriate
        % flux limiter is calculated
        %% 1-Inflow boundary
        %For dissolution, the inflow boundary is at the right of the
        %spatial domain. The ghost cells are assumed to be 0.

        f(end,n+1)=f(end,n)+dt/dL*G(n)*f(end,n)+0.5*dt/dL*G(n)*(1+dt/dL*G(n))*(fluxLimiter(end)*(-f(end,n))-fluxLimiter(end-1)*(f(end,n)-f(end-1,n)));

        %% 2-Interior volume
        f(end-1:-1:2,n+1) = f(end-1:-1:2,n) - dt/dL*(G(n)*f(end:-1:3,n)-G(n)*f(end-1:-1:2,n)) + 0.5*dt/dL*G(n)*(1+dt/dL*G(n))*(fluxLimiter(end-1:-1:2).*(f(end:-1:3,n)-f(end-1:-1:2,n))-fluxLimiter(end-2:-1:1).*(f(end-1:-1:2,n)-f(end-2:-1:1,n)));

        %% 3-Outflow boundary
        % For dissolution, the outflow boundary is at the left of the
        % spatial domain. The ghost cell is obtained using zero-order extrapolation

        % An outlet flux limiter is not required for the outflow boundary

        f(1,n+1)=f(1,n)-dt/dL*(G(n)*f(2,n)-G(n)*f(1,n))+0.5*dt/dL*G(n)*(1+dt/dL*G(n))*(fluxLimiter(1)*(f(2,n)-f(1,n)));

        %% Use liquid phase mass balance to determine supersaturation at next time step
        m3(n+1)=trapz(L,L.^3.*f(:,n+1)');
        concentration(n+1)=concentration(n)-particleDensity*shapeFactor*(m3(n+1)-m3(n));

        % Interpolate temperature to find solubility and superstaturation
        temperature(n+1) = interp1(temperatureRamp(1,:),temperatureRamp(2,:),t(n+1));
        solubility = solubilityFactor*3.37*exp(0.036*temperature(n+1));
        supersaturation(n+1)=concentration(n+1)/solubility;

        if supersaturation(n+1)>=1 % Necessary to make sure it remains a dissolution problem
            supersaturation(n+1)=1;
        end

        G(n+1)=kd1*exp(-kd2/(temperature(n+1)+273.15))*(supersaturation(n+1)-1); %[um/h];

        % Increase time counter
        n=n+1;
        
    end
end
% slice off extre cells
t = t(1:n);
m3 = m3(1:n);
concentration = concentration(1:n);
temperature = temperature(1:n);
supersaturation = supersaturation(1:n);
G = G(1:n);
f = f(:,1:n);
end