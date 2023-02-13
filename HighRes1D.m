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
%
% Purpose: Implements a high resolution finite volume method (with van Leer
% limiter) to solve for the time evolution of a particle size distribution
% over a given time range.
% References: 
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
function [f, c, G, S, m3, t, ce] = HighRes1D(dL, L, tmax, k1, k2, kv, T, ParticleDensity, c0, f0)

%Equilibrium concentration
ce = 3.37*exp(0.0359*T);

%Initial values
f(:,1)=f0;
m3(1)=trapz(L.^3.*f(:,1)');
c(1)=c0;
S(1)=c0/ce;
G(1)=k1*(S(1)-1)^k2;

%Courant number to specify the maximum stable time step
CourantNumber = 0.9;

%while loop to update and store cell averages over time
t=0;
n=1;

if G(1)>0
% Growth
    while t(n)<tmax
        
        %Calculate stable time step using available values
        dt=CourantNumber*dL/G(n);
        
        %Check if the max stable time step will exceed time range
        if tmax-t(n)<=dt
            t(n+1)=tmax;
        else
            t(n+1)=t(n)+dt;
        end
    
        %Calculate the PSD at the new time using high resolution method
        %% 1-Inflow boundary
        i = 1;
    
        %Determine smoothness and calculate appropriate flux limiter for cell
        %inlet and outlet fluxes
        
        fluxLimiterIn = 1;
        fluxLimiterOut = 1;
    
        f(i,n+1)=f(i,n)-CourantNumber*f(i,n)-0.5*CourantNumber*(1-CourantNumber)*(fluxLimiterOut*(f(i+1,n)-f(i,n))-fluxLimiterIn*(f(i,n)));
    
        %% 2-Interior volume
        for i=2:length(L)-1
    
            fluxLimiterIn = fluxLimiterOut;
            if f(i+1,n)-f(i,n)==0
                fluxLimiterOut = 1;
            else
                smoothnessOut = (f(i,n)-f(i-1,n))/(f(i+1,n)-f(i,n));
                fluxLimiterOut = (smoothnessOut+abs(smoothnessOut))/(1+abs(smoothnessOut));
            end
            
            f(i,n+1)=f(i,n)-CourantNumber*(f(i,n)-f(i-1,n))-0.5*CourantNumber*(1-CourantNumber)*(fluxLimiterOut*(f(i+1,n)-f(i,n))-fluxLimiterIn*(f(i,n)-f(i-1,n)));
        end
    
        %% 3-Outflow boundary
        i = length(L);
    
        fluxLimiterIn = fluxLimiterOut;
        % An outlet flux limiter is not required for the outflow boundary
        
        f(i,n+1)=f(i,n)-CourantNumber*(f(i,n)-f(i-1,n))-0.5*CourantNumber*(1-CourantNumber)*(fluxLimiterOut*(f(i,n)-f(i,n))-fluxLimiterIn*(f(i,n)-f(i-1,n)));
    
        % Use liquid phase mass balance to determine supersaturation at next time step 
        m3(n+1)=trapz(L.^3.*f(:,n+1)');
        c(n+1)=c(n)-ParticleDensity*kv*(m3(n+1)-m3(n));
        S(n+1)=c(n+1)/ce;
    
        if S(n+1)<=1 % Necessary to make sure it remains a growth problem
            S(n+1)=1;
        end
    
        G(n+1)=k1*(S(n+1)-1)^k2;
     
        % Increase time counter
        n=n+1;
    end
  
elseif G(1)<0
% Dissolution
% Courant number is now negative
CourantNumber = -CourantNumber;
    while t(n)<tmax
        
        %Calculate stable time step using available values
        dt=CourantNumber*dL/G(n);
        
        %Check if the max stable time step will exceed time range
        if tmax-t(n)<=dt
            t(n+1)=tmax;
        else
            t(n+1)=t(n)+dt;
        end
    
        %Calculate the PSD at the new time using high resolution method:
        % for 1, 2 & 3, the smoothness is determined and the appropriate
        % flux limiter is calculated
        %% 1-Inflow boundary
        %For dissolution, the inflow boundary is at the right of the
        %spatial domain. The ghost cells are assumed to be 0.
        i = length(L);
        
        fluxLimiterIn = 1;
        fluxLimiterOut = 1;
    
        f(i,n+1)=f(i,n)+CourantNumber*f(i,n)+0.5*CourantNumber*(1+CourantNumber)*(fluxLimiterIn*(-f(i,n))-fluxLimiterOut*(f(i,n)-f(i-1,n)));
    
        %% 2-Interior volume
        for i=length(L)-1:-1:2
    
            fluxLimiterIn = fluxLimiterOut;
            if f(i,n)-f(i-1,n)==0
                fluxLimiterOut = 1;
            else
                smoothnessOut = (f(i+1,n)-f(i,n))/(f(i,n)-f(i-1,n));
                fluxLimiterOut = (smoothnessOut+abs(smoothnessOut))/(1+abs(smoothnessOut));
            end
            
            f(i,n+1)=f(i,n)-CourantNumber*(f(i+1,n)-f(i,n))+0.5*CourantNumber*(1+CourantNumber)*(fluxLimiterIn*(f(i+1,n)-f(i,n))-fluxLimiterOut*(f(i,n)-f(i-1,n)));
        end
    
        %% 3-Outflow boundary
        % For dissolution, the outflow boundary is at the left of the
        % spatial domain. The ghost cell is obtained using zero-order extrapolation
        i = 1;
    
        fluxLimiterIn = fluxLimiterOut;
        % An outlet flux limiter is not required for the outflow boundary
        
        f(i,n+1)=f(i,n)-CourantNumber*(f(i+1,n)-f(i,n))+0.5*CourantNumber*(1+CourantNumber)*(fluxLimiterIn*(f(i+1,n)-f(i,n))-fluxLimiterOut*(f(i,n)-f(i,n)));
    
        %% Use liquid phase mass balance to determine supersaturation at next time step 
        m3(n+1)=trapz(L.^3.*f(:,n+1)');
        c(n+1)=c(n)-ParticleDensity*kv*(m3(n+1)-m3(n));
        S(n+1)=c(n+1)/ce;
    
        if S(n+1)>=1 % Necessary to make sure it remains a dissolution problem
            S(n+1)=1;
        end
    
        G(n+1)=k1*(S(n+1)-1)^k2;
     
        % Increase time counter
        n=n+1;
    end
end
end