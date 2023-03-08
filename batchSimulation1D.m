%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
% 
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
% 
% Last modified:
% - 2023/02/24, MA: initial creation
% - 2023/03/02, MA: added dialogue boxes to gather user input
% - 2023/03/06, MA: added additive dependency
% - 2023/03/07, MA: removed additive dependency (to a separate file)
% - 2023/03/08, MA: fixed trapz function
%
% Purpose: to provide a framework for different simulations of a 1D batch
% crystallization model.
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
% (4) Ramkrishna, D., 2000. Population balances : theory and applications
% to particulate systems in engineering. Academic Press.
%
% Input arguments
% 
% Output arguments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Essential parameters
clear; clc; close all

% create a dialogue box to gather essential parameters from user
prompt = {'Initial concentration [g/kg]',['k1 [' char(181) 'm]'],'k2','shape factor', ['particle density [g/' char(181) 'm' char(179) ']'],'simulation time [h]', ['length step size [' char(181) 'm]']};
dlgtitle = 'Essential parameters';
definput = {'30','10','1','1','1.46e-12','100','1'};
essentialParameters = inputdlg(prompt,dlgtitle, [1 35],definput);

initialConcentration = str2double(essentialParameters(1)); % [g/kg]

% equilibrium temperature at initial concentration (starting point of
% any simulation)
T0 = (1/0.0359)*log(initialConcentration/3.37);

% growth parameters
k1 = str2double(essentialParameters(2)); % [um/h]
k2 = str2double(essentialParameters(3));

shapeFactor = str2double(essentialParameters(4)); % assuming cuboidal particles
particleDensity = str2double(essentialParameters(5)); % [g/um3]
simulationTime = str2double(essentialParameters(6)); % [h] range of t required

% define length step and spatial domain
dL = str2double(essentialParameters(7)); % [um]
L = 1:dL:1000; % [um]

% prompt user for initial PSD
choicePSD = listdlg('ListString',{'Gaussian','Pulse'},'PromptString','Please select an inital particle size distribution','SelectionMode','single');
% use choice to generate PSD
if choicePSD == 1
    initialPSD = 1e5*normpdf(L,300,20);
elseif choicePSD == 2
    initialPSD = zeros(1,length(L));
    initialPSD(100:200)=1e3;
end

% calculate initial moments using intergation
m0i=trapz(L,initialPSD);
m1i=trapz(L,L.*initialPSD);
m2i=trapz(L,L.^2.*initialPSD);
m3i=trapz(L,L.^3.*initialPSD);

% define vector containing initial conditions for method of moments
y0 = [m0i m1i m2i m3i initialConcentration];

% prompt user for operation mode
operationMode = listdlg("ListString",{'constant temperature (crash cooling)','constant supersaturation'},"PromptString",'Please select the required simulation',"SelectionMode","single");
switch operationMode
%% Constant temperature comparison

    case 1
        % create a dialogue box to get constant temperature from user
        prompt = {['Please input the constant temperature at which crystallisation takes place [' char(176) 'C]']};
        dlgtitle = 'Constant temperature';
        definput = {'25'};
        constantTemperature = str2double(inputdlg(prompt,dlgtitle, [1 35],definput));
        
        % create crash cooling temperature ramp and generate a simulation starting
        % from equilibrium then crash cooling
        temperatureRamp = [0 5 5.00001 simulationTime; T0 T0 constantTemperature ...
            constantTemperature];
        [f, concentration, G, supersaturation, m3, t, temperature] ...
            = highRes1D(dL, L, simulationTime, k1, k2, shapeFactor, ...
            temperatureRamp, particleDensity, initialConcentration, ...
            initialPSD, operationMode);
        %calculate average length for comparison to method of moments
        m0 = trapz(L,f);
        m1 = trapz(L,L'.*f);
        averageLength = m1./m0;
        
        % use less realistic temperature ramp generate an alternate simulation at a
        % constant temperature (old method) 
        temperatureRamp_alt = [0 simulationTime; constantTemperature constantTemperature];
        [f_alt, concentration_alt, G_alt, supersaturation_alt, m3_alt, t_alt, temperature_alt] ...
            = highRes1D(dL, L, simulationTime, k1, k2, shapeFactor, ...
            temperatureRamp_alt, particleDensity, initialConcentration, ...
            initialPSD, operationMode);
    
        %calculate average length for comparison to method of moments
        m0_alt = trapz(L,f_alt);
        m1_alt = trapz(L,L'.*f_alt);
        averageLength_alt = m1_alt./m0_alt;
        
        % solve using method of moments to check concentration and moments from
        % finite volume
           
        % integrate method of moments ODE system
        [t_mom, y] = ode23s(@(t_mom, y)mom(t_mom, y, temperatureRamp, k1, k2, ...
            particleDensity, shapeFactor,operationMode), [0 simulationTime], y0);
        
        % extract concentration profile predicted by method of moments
        concentration_mom = y(:,5);
        
        % extract temperature profile predicted by method of moments
        temperature_mom = interp1(temperatureRamp(1,:),temperatureRamp(2,:),t_mom);
    
        % average length predicted by method of moments
        averageLength_mom  = y(:,2)./y(:,1);
        
        % solubility curve plots
        figure(1)
        % plot of system starting at equilibrium then crash cooling
        plot(temperature,concentration,'linewidth',1.6)
        hold on
        
        % plot of system operating at constant temperature
        plot(temperature_alt, concentration_alt,':','linewidth',1.2)
        
        % plot predicted by method of moments
         plot(temperature_mom,concentration_mom,'--','linewidth',1.2)
        
        % generate solubility curve for reference
        T = 0:0.1:70;
        solubility = 3.37*exp(0.0359*T);
        
        % plot of solubility curve
        plot(T,solubility,'linewidth',1.2)
        
        set(gca,'FontSize',18)
        title('Solubility curve')
        xlabel({'Temperature', ['[' char(176) 'C]']})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution (starting from equilibrium)', ...
            'High resolution (constant temperature)', ...
            'Method of moments', ...
            'Solubility curve (equilibrium)')
        
        % concentration and average length profiles
        figure(2)
        
        % concentration profile
        subplot(2,1,1)
        
        % concentration profile of system starting at equilibrium then crash
        % cooling
        plot(t,concentration,'linewidth',1.2)
        hold on
        
        % concentration profile of system operating at constant temperature
         plot(t_alt,concentration_alt,':','linewidth',1.2)
        
        % concentration profile predicted by method of moments
        plot(t_mom,concentration_mom,'--','linewidth',1.2)
        
        set(gca,'FontSize',18)
        title('Concentration profile')
        xlabel({'time' '[h]'})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution (starting from equilibrium)', ...
            'High resolution (constant temperature)', ...
            'Method of moments')
        
        % average length profile
        subplot(2,1,2)
        % average length profile of system starting at equilibrium then crash
        % cooling:
        plot(t,averageLength,'linewidth',1.2)
        hold on
        % average length profile of system operating at constant temperature
        plot(t_alt,averageLength_alt,':','linewidth',1.2)
        % average length profile predicted by method of moments
        plot(t_mom,averageLength_mom,'--','linewidth',1.2)
    
        set(gca,'FontSize',18)  
        title('Average length profile')
        xlabel({'time' '[h]'})
        ylabel({'Average length', ['[' char(181) 'm]']})
        legend('High resolution (starting from equilibrium)', ...
            'High resolution (constant temperature)', ...
            'Method of moments')
    
        % PSD plots:
        figure(3)
        % PSD evolution comparison
        plot(L,f(:,1),'linewidth',1.2)
        hold on
        plot(L,f(:,end),'linewidth',1.2)
    
        plot(L,f_alt(:,1),'--','linewidth',1.2)
        hold on
        plot(L,f_alt(:,end),'--','linewidth',1.2)
        
        set(gca,'FontSize',18)  
        title('PSD evolution comparison')
        xlabel({'Length', ['[' char(181) 'm]']})
        ylabel({'Number density', '[μm^{-1} kg^{-1}]'})
        legend('Initial PSD (starting from equilibrium)', ...
            'Final PSD (starting from equilibrium)', ...
            'Initial PSD (constant temperature)', ...
            'Final PSD (constant temperature')

%% Constant supersaturation
    case 2
        % create dialogue box to get constant supersaturation from user
        prompt = {'Constant supersaturation','Final temperature (equilibrium)'};
        dlgtitle = 'Constant supersaturation';
        definput = {'1.5','25'};
        answer = inputdlg(prompt,dlgtitle, [1 35],definput);
        constantSupersaturation = str2double(answer(1));
        finalTemperature = str2double(answer(2));
    
        % generate temperature ramp by calculating temperature at each
        % concentration that ensures a constant supersaturation
        finalconcentration = 3.37*constantSupersaturation*exp(0.0359*finalTemperature);
        if finalconcentration > initialConcentration
            % in the case of dissolution
            c = initialConcentration:0.1:finalconcentration;
            T = (1/0.0359)*log(c/(3.37*constantSupersaturation));
            temperatureRamp = [c, finalconcentration+5; T ,finalTemperature]; % dummy cell added to prevent concentrations higher than final concentration causing errors
        else
            % in the case of growth
            c = finalconcentration:0.1:initialConcentration+1; 
            T = (1/0.0359)*log(c/(3.37*constantSupersaturation));
            temperatureRamp = [finalconcentration-10, c;finalTemperature, T]; % dummy cell added to prevent concentrations lower than final concentration causing errors
        end
    
        % High resolution solution:
        [f, concentration, G, supersaturation, m3, t, temperature] ...
            = highRes1D(dL, L, simulationTime, k1, k2, shapeFactor, ...
            temperatureRamp, particleDensity, initialConcentration, ...
            initialPSD, operationMode);
        %calculate average length for comparison to method of moments
        m0 = trapz(L,f);
        m1 = trapz(L,L'.*f);
        averageLength = m1./m0;
        
        % method of moments solution:
        [t_mom, y] = ode23s(@(t_mom, y)mom(t_mom, y, temperatureRamp, k1, k2, ...
            particleDensity, shapeFactor,operationMode), [0 simulationTime], y0);
        % extract concentration profile predicted by method of moments
        concentration_mom = y(:,5);
        % average length predicted by method of moments
        averageLength_mom  = y(:,2)./y(:,1);
        
        % extract temperature profile predicted by method of moments
        temperature_mom = interp1(temperatureRamp(1,:),temperatureRamp(2,:),concentration_mom);
    
        % solubility curve plots
        figure(1)
        % plot of system starting at equilibrium then crash cooling
        plot(temperature,concentration,'linewidth',1.6)
        hold on
        
        % plot predicted by method of moments
    
         plot([T0; temperature_mom],[initialConcentration; concentration_mom],'--','linewidth',1.2)
        
        % generate solubility curve for reference
        T = 0:0.1:70;
        solubility = 3.37*exp(0.0359*T);
        
        % plot of solubility curve to show driving force
        plot(T,solubility,'linewidth',1.2)
        
        set(gca,'FontSize',18)
        title('Solubility curve')
        xlabel({'Temperature', ['[' char(176) 'C]']})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution', ...
            'Method of moments ', ...
            'Solubility curve')
        
        % concentration and average length plots
        figure(2)
        subplot(2,1,1)
        
        % concentration profile of system starting at equilibrium then crash
        % cooling
        plot(t,concentration,'linewidth',1.2)
        hold on
        
        % concentration profile predicted by method of moments
        plot(t_mom,concentration_mom,'--','linewidth',1.2)
        
        set(gca,'FontSize',18)
        title('Concentration profile')
        xlabel({'time' '[h]'})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution','Method of moments')
        
        % average length profile
        subplot(2,1,2)
        % average length profile of system starting at equilibrium then crash
        % cooling:
        plot(t,averageLength,'linewidth',1.2)
        hold on
        % average length profile predicted by method of moments:
        plot(t_mom,averageLength_mom,'--','linewidth',1.2)
    
        set(gca,'FontSize',18)  
        title('Average length profile')
        xlabel({'time' '[h]'})
        ylabel({'Average length', ['[' char(181) ']']})
        legend('High resolution','Method of moments')

         % PSD plots:
        figure(3)
        % PSD evolution comparison
        plot(L,f(:,1),'linewidth',1.2)
        hold on
        plot(L,f(:,end),'linewidth',1.2)
        
        set(gca,'FontSize',18)  
        title('PSD evolution comparison')
        xlabel({'time' '[h]'})
        ylabel({'Number density', '[μm^{-1} kg^{-1}]'})
        legend('Initial PSD (starting from equilibrium)', ...
            'Final PSD (starting from equilibrium)')
    % PSD plots:

    %% Custom temperature ramp
% elseif operationMode == 3
% prompt user for temperature ramp
% temperatureRamp = input(['\nProvide the desired temperature ramp as a matrix' ...
%     ' over the entire timespan where the top row represents time while' ...
%     'the bottom row represents corresponding temperature.\n(for a constant' ...
%     ' temperature, T, simply type [0 100; T T]\n']);
end