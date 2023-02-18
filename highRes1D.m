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
% - 2023/02/16, MA: Added temperature dependence
% - 2023/02/16, MA: Replaced spatial for loop with vector operations
%
% Purpose: Implements a high resolution finite volume method (with van Leer
% limiter) to solve for the time evolution of a particle size distribution
% over a given time range.
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
% kv: Scalar representing particle shape factor
%
% T: Scalar representing the temperature
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
function [f, concentration, G, supersaturation, m3, t, solubility] = highRes1D(dL, L, simulationTime, k1, k2, shapeFactor, temperatureRamp, ParticleDensity, initialConcentration, f0)

%Equilibrium concentration


%Initial values
f(:,1)=f0;
m3(1)=trapz(L.^3.*f(:,1)');
concentration(1)=initialConcentration;
solubility = 3.37*exp(0.0359*temperatureRamp(2,1));
supersaturation(1)=initialConcentration/solubility;
G(1)=k1*(supersaturation(1)-1)^k2;

%initialise smoothness and flux limiter functions
smoothness = zeros(length(L),1);
fluxLimiter = zeros(length(L),1);

%Courant number to specify the maximum stable time step
CourantNumber = 0.9;

%while loop to update and store cell averages over time
t=0;
n=1;

if G(1)>0
% Growth
    while t(n)<simulationTime
        
        %Calculate stable time step using available values
        dt=CourantNumber*dL/G(n);
        
        %Check if the max stable time step will exceed time range
        if simulationTime-t(n)<=dt
            t(n+1)=simulationTime;
        else
            t(n+1)=t(n)+dt;
        end
    
        %Calculate the PSD at the new time using high resolution method
        %% 1-Inflow boundary
    
        %Determine smoothness and calculate appropriate flux limiter for cell
        %inlet and outlet fluxes
        
        fluxLimiter(1) = 1;
        fluxLimiter(2) = 1;
    
        f(1,n+1)=f(1,n)-CourantNumber*f(1,n)-0.5*CourantNumber*(1-CourantNumber)*(fluxLimiter(2)*(f(2,n)-f(1,n))-fluxLimiter(1)*(f(1,n)));
    
        %% 2-Interior volume

        smoothness(3:end) = (f(2:end-1,n)-f(1:end-2,n)+eps)./(f(3:end,n)-f(2:end-1,n)+eps);
        fluxLimiter(3:end) = (smoothness(3:end)+abs(smoothness(3:end)))./(1+abs(smoothness(3:end)));
        f(2:end-1,n+1) = f(2:end-1,n) - CourantNumber*(f(2:end-1,n)-f(1:end-2,n)) - 0.5*CourantNumber*(1-CourantNumber)*(fluxLimiter(3:end).*(f(3:end,n)-f(2:end-1,n)) - fluxLimiter(2:end-1).*(f(2:end-1,n)-f(1:end-2,n)));

        %% 3-Outflow boundary

        % An outlet flux limiter is not required for the outflow boundary
        
        f(end,n+1)=f(end,n)-CourantNumber*(f(end,n)-f(end-1,n))-0.5*CourantNumber*(1-CourantNumber)*(-fluxLimiter(end)*(f(end,n)-f(end-1,n)));
    
        %% Use liquid phase mass balance to determine supersaturation at next time step 
        m3(n+1)=trapz(L.^3.*f(:,n+1)');
        concentration(n+1)=concentration(n)-ParticleDensity*shapeFactor*(m3(n+1)-m3(n));
        
        % Interpolate temperature to find solubility ans superstaturation
        solubility = 3.37*exp(0.0359*interp1(temperatureRamp(1,:),temperatureRamp(2,:),t(n+1)));
        supersaturation(n+1)=concentration(n+1)/solubility;
    
        if supersaturation(n+1)<=1 % Necessary to make sure it remains a growth problem
            supersaturation(n+1)=1;
        end
    
        G(n+1)=k1*(supersaturation(n+1)-1)^k2;
     
        % Increase time counter
        n=n+1;
    end
  
elseif G(1)<0
% Dissolution
% Courant number is now negative
CourantNumber = -CourantNumber;
    while t(n)<simulationTime
        
        %Calculate stable time step using available values
        dt=CourantNumber*dL/G(n);
        
        %Check if the max stable time step will exceed time range
        if simulationTime-t(n)<=dt
            t(n+1)=simulationTime;
        else
            t(n+1)=t(n)+dt;
        end
    
        %Calculate the PSD at the new time using high resolution method:
        % for 1, 2 & 3, the smoothness is determined and the appropriate
        % flux limiter is calculated
        %% 1-Inflow boundary
        %For dissolution, the inflow boundary is at the right of the
        %spatial domain. The ghost cells are assumed to be 0.
        
        fluxLimiter(end) = 1;
        fluxLimiter(end-1) = 1;
    
        f(end,n+1)=f(end,n)+CourantNumber*f(end,n)+0.5*CourantNumber*(1+CourantNumber)*(fluxLimiter(end)*(-f(end,n))-fluxLimiter(end-1)*(f(end,n)-f(end-1,n)));
    
        %% 2-Interior volume
%         for i=length(L)-1:-1:2
%     
%             fluxLimiterIn = fluxLimiterOut;
%             if f(i,n)-f(i-1,n)==0
%                 fluxLimiterOut = 1;
%             else
%                 smoothnessOut = (f(i+1,n)-f(i,n))/(f(i,n)-f(i-1,n));
%                 fluxLimiterOut = (smoothnessOut+abs(smoothnessOut))/(1+abs(smoothnessOut));
%             end
%             
%             f(i,n+1)=f(i,n)-CourantNumber*(f(i+1,n)-f(i,n))+0.5*CourantNumber*(1+CourantNumber)*(fluxLimiterIn*(f(i+1,n)-f(i,n))-fluxLimiterOut*(f(i,n)-f(i-1,n)));
%         end
    
        smoothness(end-2:-1:1) = (f(end:-1:3,n)-f(end-1:-1:2,n)+eps)./(f(end-1:-1:2,n)-f(end-2:-1:1,n)+eps);
        fluxLimiter(end-2:-1:1) = (smoothness(end-2:-1:1)+abs(smoothness(end-2:-1:1)))./(1+abs(smoothness(end-2:-1:1)));
        f(end-1:-1:2,n+1) = f(end-1:-1:2,n) - CourantNumber*(f(end:-1:3,n)-f(end-1:-1:2,n)) + 0.5*CourantNumber*(1+CourantNumber)*(fluxLimiter(end-2:-1:1).*(f(end:-1:3,n)-f(end-1:-1:2,n))-fluxLimiter(end-1:-1:2).*(f(end-1:-1:2,n)-f(end-2:-1:1,n)));

        %% 3-Outflow boundary
        % For dissolution, the outflow boundary is at the left of the
        % spatial domain. The ghost cell is obtained using zero-order extrapolation

        % An outlet flux limiter is not required for the outflow boundary
        
        f(1,n+1)=f(1,n)-CourantNumber*(f(2,n)-f(1,n))+0.5*CourantNumber*(1+CourantNumber)*(fluxLimiter(1)*(f(1+1,n)-f(1,n)));
    
        %% Use liquid phase mass balance to determine supersaturation at next time step 
        m3(n+1)=trapz(L.^3.*f(:,n+1)');
        concentration(n+1)=concentration(n)-ParticleDensity*shapeFactor*(m3(n+1)-m3(n));

         % Interpolate temperature to find solubility ans superstaturation
        solubility = 3.37*exp(0.0359*interp1(temperatureRamp(1,:),temperatureRamp(2,:),t(n+1)));
        supersaturation(n+1)=concentration(n+1)/solubility;
    
        if supersaturation(n+1)>=1 % Necessary to make sure it remains a dissolution problem
            supersaturation(n+1)=1;
        end
    
        G(n+1)=k1*(supersaturation(n+1)-1)^k2;
     
        % Increase time counter
        n=n+1;
    end
end
end