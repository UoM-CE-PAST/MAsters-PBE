%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
%
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
%
% Last modified:
% - 2023/03/22, MA: initial creation
%
%
% Purpose: Implements a high resolution finite volume method (with van Leer
% limiter) to solve for the time evolution of a 2D particle size distribution
% over a given time range. This function differentiates itself by using the
% solubility and growth rate outlined in the paper by T. Vetter et al. This
% version is designed to be a better performance version of
% highRes2D_additive.m by only storing the initial and final PSSDs.
%
% References:
% (1) LeVeque, R.J., 2002. Finite Volume Methods for Hyperbolic Problems, Cambridge Texts in Applied Mathematics. Cambridge University Press, Cambridge. https://doi.org/10.1017/CBO9780511791253
% (2) Gunawan, R., Fusman, I., Braatz, R.D., 2004. High resolution algorithms for multidimensional population balance equations. AIChE Journal 50. https://doi.org/10.1002/aic.10228
% (3) Ma, D.L., Tafti, D.K., Braatz, R.D., 2002. High-resolution simulation of multidimensional crystal growth. Industrial and Engineering Chemistry Research 41. https://doi.org/10.1021/ie010680u
% (4) Vetter, T., Mazzotti, M., Brozio, J., 2011. Slowing the growth rate
% of ibuprofen crystals using the polymeric additive pluronic F127. Crystal
% Growth and Design 11. https://doi.org/10.1021/cg200352u
%
% Input arguments
%
% Output arguments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveTD(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, temperatureRamp, particleDensity, initialConcentration, initialPSD,solubilityFactor,growthFactor)

%initialise variables (overspecify then slice off later)
t = zeros(1,50000);
m00 = zeros(1,length(t));
m31 = zeros(1,length(t));
m22 = zeros(1,length(t));
m21 = zeros(1,length(t));
m41 = zeros(1,length(t));
m23 = zeros(1,length(t));

concentration = zeros(1,length(t));
temperature = zeros(1,length(t));
supersaturation = zeros(1,length(t));
G1 = zeros(1,length(t));
G2 = zeros(1,length(t));
f_mid = zeros(length(L2),length(L1)); %initialise intermediate solution

% Initial values
f_old=initialPSD;
m00(1) = sum(initialPSD,'all')*dL1*dL2;
m31(1) = sum(L1.^3.*L2'.*initialPSD,'all')*dL1*dL2;
m22(1) = sum(L1.^2.*L2'.^2.*initialPSD,'all')*dL1*dL2;
m21(1) = sum(L1.^2.*L2'.*initialPSD,'all')*dL1*dL2;
m41(1) = sum(L1.^4.*L2'.*initialPSD,'all')*dL1*dL2;
m23(1) = sum(L1.^2.*L2'.^3.*initialPSD,'all')*dL1*dL2;
concentration(1) = initialConcentration;
T0 = (1/0.036)*log(initialConcentration/(3.37*solubilityFactor));
temperature(1) = T0;
solubility = solubilityFactor*3.37*exp(0.036*temperature(1));
supersaturation(1)=initialConcentration/solubility;
if supersaturation(1) > 1
    G1(1) = kg11*exp(-kg12/(temperature(1)+273.15))*(supersaturation(1)-1)^kg13;
    G2(1) = growthFactor*kg21*exp(-kg22/(temperature(1)+273.15))*(supersaturation(1)-1)^kg23;
elseif supersaturation(1) < 1
    G1(1) = kd11*exp(-kd12/(temperature(1)+273.15))*(supersaturation(1)-1);
    G2(1) = kd21*exp(-kd22/(temperature(1)+273.15))*(supersaturation(1)-1);
else
    G1(1) = 0;
    G2(2) = 0;
end

% initialise smoothness and flux limiter functions
smoothness = zeros(length(L2),length(L1));
fluxLimiter = zeros(length(L2),length(L1));

%Courant number to specify the maximum stable time step
CourantNumber = 1;

% set time index
n=1;

% if system starts at equilibrium (i.e. G=0) then ensure system doesn't
% change until driving force is non-zero
while temperature(n) == T0 || supersaturation(n) == 1
    t(n+1)=temperatureRamp(1,n+1);
    f = initialPSD;
    m00(n+1) = sum(initialPSD,'all')*dL1*dL2;
    m31(n+1) = sum(L1.^3.*L2'.*initialPSD,'all')*dL1*dL2;
    m22(n+1) = sum(L1.^2.*L2'.^2.*initialPSD,'all')*dL1*dL2;
    m21(n+1) = sum(L1.^2.*L2'.*initialPSD,'all')*dL1*dL2;
    m41(n+1) = sum(L1.^4.*L2'.*initialPSD,'all')*dL1*dL2;
    m23(n+1) = sum(L1.^2.*L2'.^3.*initialPSD,'all')*dL1*dL2;
    concentration(n+1)=initialConcentration;
    temperature(n+1) = temperatureRamp(2,n+1);
    solubility = solubilityFactor*3.37*exp(0.036*temperature(n+1));
    supersaturation(n+1)=concentration(n+1)/solubility;supersaturation(n+1)=initialConcentration/solubility;

    if supersaturation(n+1)>1
        G1(n+1) = kg11*exp(-kg12/(temperature(n+1)+273.15))*(supersaturation(n+1)-1)^kg13;
        G2(n+1) = growthFactor*kg21*exp(-kg22/(temperature(n+1)+273.15))*(supersaturation(n+1)-1)^kg23;
    elseif supersaturation(n+1)<1
        G1(n+1) = kd11*exp(-kd12/(temperature(n+1)+273.15))*(supersaturation(n+1)-1);
        G2(n+1) = kd21*exp(-kd22/(temperature(n+1)+273.15))*(supersaturation(n+1)-1);
    else
        G1(n+1) = 0;
        G2(n+1) = 0;
    end
    n = n+1;
    f_old = f;
end

%while loop to update and store cell averages over time
while t(n)<simulationTime
    if supersaturation(n) >= 1
        % Growth
        % Calculate stable time step using available values
        dt=CourantNumber/(G1(n)/dL1 + G2(n)/dL2);

        %Check if the max stable time step will exceed time range
        if simulationTime-t(n)<=dt
            t(n+1)=simulationTime;
            dt = simulationTime-t(n);
        else
            t(n+1)=t(n)+dt;
        end

        %% L1 sweeps
        j = 1:length(L2);
        % Update the PSD at the new time using high resolution method
        %% 1-Inflow boundary

        %Determine smoothness and calculate appropriate flux limiter for cell
        %inlet and outlet fluxes

        fluxLimiter(j,1) = 1;
        fluxLimiter(j,2) = 1;

        f_mid(j,1)=f_old(j,1)-dt/dL1*G1(n)*f_old(j,1)-0.5*dt/dL1*G1(n)*(1-dt/dL1*G1(n))*(fluxLimiter(j,2).*(f_old(j,2)-f_old(j,1))-fluxLimiter(j,1).*(f_old(j,1)));

        %% 2-Interior volume

        smoothness(j,3:end) = (f_old(j,2:end-1)-f_old(j,1:end-2)+eps)./(f_old(j,3:end)-f_old(j,2:end-1)+eps);
        fluxLimiter(j,3:end) = (smoothness(j,3:end)+abs(smoothness(j,3:end)))./(1+abs(smoothness(j,3:end)));
        f_mid(j,2:end-1) = f_old(j,2:end-1) - dt/dL1*(G1(n)*f_old(j,2:end-1)-G1(n)*f_old(j,1:end-2)) - 0.5*dt/dL1*G1(n)*(1-dt/dL1*G1(n))*(fluxLimiter(j,3:end).*(f_old(j,3:end)-f_old(j,2:end-1)) - fluxLimiter(j,2:end-1).*(f_old(j,2:end-1)-f_old(j,1:end-2)));

        %% 3-Outflow boundary

        % An outlet flux limiter is not required for the outflow boundary

        f_mid(j,end)=f_old(j,end)-dt/dL1*(G1(n)*f_old(j,end)-G1(n)*f_old(j,end-1))-0.5*dt/dL1*G1(n)*(1-dt/dL1*G1(n))*(-fluxLimiter(j,end).*(f_old(j,end)-f_old(j,end-1)));

        %% L2 sweeps
        i = 1:length(L1);
        %% 1-Inflow boundary

        %Determine smoothness and calculate appropriate flux limiter for cell
        %inlet and outlet fluxes

        fluxLimiter(1,i) = 1;
        fluxLimiter(2,i) = 1;

        f(1,i)=f_mid(1,i)-dt/dL2*G2(n)*f_mid(1,i)-0.5*dt/dL2*G2(n)*(1-dt/dL2*G2(n))*(fluxLimiter(2,i).*(f_mid(2,i)-f_mid(1,i))-fluxLimiter(1,i).*(f_mid(1,i)));

        %% 2-Interior volume

        smoothness(3:end,i) = (f_mid(2:end-1,i)-f_mid(1:end-2,i)+eps)./(f_mid(3:end,i)-f_mid(2:end-1,i)+eps);
        fluxLimiter(3:end,i) = (smoothness(3:end,i)+abs(smoothness(3:end,i)))./(1+abs(smoothness(3:end,i)));
        f(2:end-1,i) = f_mid(2:end-1,i) - dt/dL2*(G2(n)*f_mid(2:end-1,i)-G2(n)*f_mid(1:end-2,i)) - 0.5*dt/dL2*G2(n)*(1-dt/dL2*G2(n))*(fluxLimiter(3:end,i).*(f_mid(3:end,i)-f_mid(2:end-1,i)) - fluxLimiter(2:end-1,i).*(f_mid(2:end-1,i)-f_mid(1:end-2,i)));

        %% 3-Outflow boundary

        % An outlet flux limiter is not required for the outflow boundary

        f(end,i)=f_mid(end,i)-dt/dL2*(G2(n)*f_mid(end,i)-G2(n)*f_mid(end-1,i))-0.5*dt/dL2*G2(n)*(1-dt/dL2*G2(n))*(-fluxLimiter(end,i).*(f_mid(end,i)-f_mid(end-1,i)));


        %% Use liquid phase mass balance to determine supersaturation at next time step
        m00(n+1) = sum(f,'all')*dL1*dL2;
        m31(n+1) = sum(L1.^3.*L2'.*f,'all')*dL1*dL2;
        m22(n+1) = sum(L1.^2.*L2'.^2.*f,'all')*dL1*dL2;
        m21(n+1)=sum(L1.^2.*L2'.*f,'all')*dL1*dL2;
        m41(n+1) = sum(L1.^4.*L2'.*f,'all')*dL1*dL2;
        m23(n+1) = sum(L1.^2.*L2'.^3.*f,'all')*dL1*dL2;
        concentration(n+1)=concentration(n)-particleDensity*shapeFactor*(m21(n+1)-m21(n));

        % Interpolate temperature to find solubility and superstaturation

        % temperature is time-dependent, interpolate with time
        temperature(n+1) = interp1(temperatureRamp(1,:),temperatureRamp(2,:),t(n+1));
        solubility = solubilityFactor*3.37*exp(0.036*temperature(n+1));
        supersaturation(n+1)=concentration(n+1)/solubility;

        if supersaturation(n+1) >= 1
            G1(n+1) = kg11*exp(-kg12/(temperature(n+1)+273.15))*(supersaturation(n+1)-1)^kg13;
            G2(n+1) = growthFactor*kg21*exp(-kg22/(temperature(n+1)+273.15))*(supersaturation(n+1)-1)^kg23;
        elseif supersaturation(n+1) < 1
            G1(n+1) = kd11*exp(-kd12/(temperature(n+1)+273.15))*(supersaturation(n+1)-1);
            G2(n+1) = kd21*exp(-kd22/(temperature(n+1)+273.15))*(supersaturation(n+1)-1);
        end
        f_old = f;
        % Increase time counter
        n=n+1;
    elseif supersaturation(n) <= 1
        % Dissolution
        % Courant number is now negative
        dt=-CourantNumber/(G1(n)/dL1 + G2(n)/dL2);

        %Check if the max stable time step will exceed time range
        if simulationTime-t(n)<=dt
            t(n+1)=simulationTime;
            dt = simulationTime-t(n);
        else
            t(n+1)=t(n)+dt;
        end

        % Calculate the PSD at the new time using the high resolution
        % method: for 1, 2 & 3, the smoothness is determined and the
        % appropriate flux limiter is calculated

        %% L1 sweeps
        j = 1:length(L2);
        %% 1-Inflow boundary
        % For dissolution, the inflow boundary is at the right of the
        % spatial domain. The ghost celss are assumed to be 0.

        fluxLimiter(j,end) = 1;
        fluxLimiter(j,end) = 1;

        f_mid(j,end)=f_old(j,end)+dt/dL1*G1(n)*f_old(j,end)+0.5*dt/dL1*G1(n)*(1+dt/dL1*G1(n))*(fluxLimiter(j,end).*(-f_old(j,end))-fluxLimiter(j,end-1).*(f_old(j,end)-f_old(end-1)));

        %% 2-Interior volume
        smoothness(j,end-2:-1:1) = (f_old(j,end:-1:3)-f_old(j,end-1:-1:2)+eps)./(f_old(j,end-1:-1:2)-f_old(j,end-2:-1:1)+eps);
        fluxLimiter(j,end-2:-1:1) = (smoothness(j,end-2:-1:1)+abs(smoothness(j,end-2:-1:1)))./(1+abs(smoothness(j,end-2:-1:1)));
        f_mid(j,end-1:-1:2) = f_old(j,end-1:-1:2) - dt/dL1*(G1(n)*f_old(j,end:-1:3)-G1(n)*f_old(j,end-1:-1:2)) + 0.5*dt/dL1*G1(n)*(1+dt/dL1*G1(n))*(fluxLimiter(j,end-2:-1:1).*(f_old(j,end:-1:3)-f_old(j,end-1:-1:2)) - fluxLimiter(j,end-1:-1:2).*(f_old(j,end-1:-1:2)-f_old(j,end-2:-1:1)));

        %% 3-Outflow boundary
        % For dissolution, the outflow boundary is at the left of the
        % spatial domain. The ghost cell is obtained using zero-order extrapolation

        % An outlet flux limiter is not required for the outflow boundary
        f_mid(j,1)=f_old(j,1)-dt/dL1*(G1(n)*f_old(j,2)-G1(n)*f_old(j,1))+0.5*dt/dL1*G1(n)*(1+dt/dL1*G1(n))*(fluxLimiter(j,1).*(f_old(j,2)-f_old(j,1)));
        %% L2 sweeps
        i = 1:length(L1);
        %% 1-Inflow boundary

        %Determine smoothness and calculate appropriate flux limiter for cell
        %inlet and outlet fluxes

        fluxLimiter(end,i) = 1;
        fluxLimiter(end-1,i) = 1;

        f(end,i)=f_mid(end,i)+dt/dL2*G2(n)*f_mid(end,i)+0.5*dt/dL2*G2(n)*(1+dt/dL2*G2(n))*(fluxLimiter(end,i).*(-f_mid(end,i))-fluxLimiter(end-1,i).*(f_mid(end,i)-f_mid(end-1,i)));

        %% 2-Interior volume

        smoothness(end-2:-1:1,i) = (f_mid(end:-1:3,i)-f_mid(end-1:-1:2,i)+eps)./(f_mid(end-1:-1:2,i)-f_mid(end-2:-1:1,i)+eps);
        fluxLimiter(end-2:-1:1,i) = (smoothness(end-2:-1:1,i)+abs(smoothness(end-2:-1:1,i)))./(1+abs(smoothness(end-2:-1:1,i)));
        f(end-1:-1:2,i) = f_mid(end-1:-1:2,i) - dt/dL2*(G2(n)*f_mid(end:-1:3,i)-G2(n)*f_mid(end-1:-1:2,i)) + 0.5*dt/dL2*G2(n)*(1+dt/dL2*G2(n))*(fluxLimiter(end-2:-1:1,i).*(f_mid(end:-1:3,i)-f_mid(end-1:-1:2,i)) - fluxLimiter(end-1:-1:2,i).*(f_mid(end-1:-1:2,i)-f_mid(end-2:-1:1,i)));

        %% 3-Outflow boundary

        % An outlet flux limiter is not required for the outflow boundary

        f(1,i)=f_mid(1,i)-dt/dL2*(G2(n)*f_mid(2,i)-G2(n)*f_mid(1,i))+0.5*dt/dL2*G2(n)*(1+dt/dL2*G2(n))*(fluxLimiter(1,i).*(f_mid(2,i)-f_mid(1,i)));

        %% Use liquid phase mass balance to determine supersaturation at next time step
        m00(n+1) = sum(f,'all')*dL1*dL2;
        m31(n+1) = sum(L1.^3.*L2'.*f,'all')*dL1*dL2;
        m22(n+1) = sum(L1.^2.*L2'.^2.*f,'all')*dL1*dL2;
        m21(n+1)=sum(L1.^2.*L2'.*f,'all')*dL1*dL2;
        m41(n+1) = sum(L1.^4.*L2'.*f,'all')*dL1*dL2;
        m23(n+1) = sum(L1.^2.*L2'.^3.*f,'all')*dL1*dL2;
        concentration(n+1)=concentration(n)-particleDensity*shapeFactor*(m21(n+1)-m21(n));
        % temperature is time-dependent, interpolate with time
        temperature(n+1) = interp1(temperatureRamp(1,:),temperatureRamp(2,:), t(n+1));
        solubility = solubilityFactor*3.37*exp(0.036*temperature(n+1));
        supersaturation(n+1)=concentration(n+1)/solubility;

        if supersaturation(n+1) <= 1
            G1(n+1) = kd11*exp(-kd12/(temperature(n+1)+273.15))*(supersaturation(n+1)-1);
            G2(n+1) = kd21*exp(-kd22/(temperature(n+1)+273.15))*(supersaturation(n+1)-1);
        elseif supersaturation(n+1) > 1
            G1(n+1) = kg11*exp(-kg12/(temperature(n+1)+273.15))*(supersaturation(n+1)-1)^kg13;
            G2(n+1) = growthFactor*kg21*exp(-kg22/(temperature(n+1)+273.15))*(supersaturation(n+1)-1)^kg23;
        end
        f_old = f;
        n=n+1;
        % Increase time counter

    end
end
t = t(1:n);
m00 = m00(1:n);
m31 = m31(1:n);
m22 = m22(1:n);
m21 = m21(1:n);
m41 = m41(1:n);
m23 = m23(1:n);
concentration = concentration(1:n);
temperature = temperature(1:n);
supersaturation = supersaturation(1:n);
G1 = G1(1:n);
G2 = G2(1:n);
end