%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The University of Manchester, United Kingdom
% 
% Project: MEng Dissertation
% Year: 2023
% MATLAB: R2022b
% Author(s): Mohammed Alsubeihi (MA)
% 
% Last modified:
% - 2023/03/07, MA: initial creation
% - 2023/03/08, MA: fixed trapz function
%
% Purpose: to provide a framework for different simulations of a 1D batch
% crystallization model in the presence of additive.
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
% (5) Vetter, T., Mazzotti, M., Brozio, J., 2011. Slowing the growth rate
% of ibuprofen crystals using the polymeric additive pluronic F127. Crystal
% Growth and Design 11. https://doi.org/10.1021/cg200352u
%
% Input arguments
% 
% Output arguments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Essential parameters
clear; clc; close all

% create a dialogue box to gather essential parameters from user
prompt = {'Initial concentration [g/kg]','additive concentration [kg/kg]s','shape factor', ['particle density [g/' char(181) 'm' char(179) ']'],'simulation time [h]', ['length step size [' char(181) 'm]']};
dlgtitle = 'Essential parameters';
definput = {'50','0','0.5236','1.1128e-12','9000','1'};
essentialParameters = inputdlg(prompt,dlgtitle, [1 35],definput);

initialConcentration = str2double(essentialParameters(1)); % [g/kg]

additiveConcentration = str2double(essentialParameters(2)); % [kg/kg]

% array containing the solubility parameters for each additive
% concentration: (c; p4; p3; p2; p1; p0) 
solubilityParameters = [0 0.04 0.08; 3.594e-3 2.719e-3 2.408e-3;
    9.957e-3 8.747e-3 9.19e-3;
    1.377e-2 1.785e-2 2.039e-2;
    2.841e-2 3.658e-2 4.265e-2;
    4.171e-2 5.046e-2 6.055e-2];
p4 = interp1(solubilityParameters(1,:),solubilityParameters(2,:),additiveConcentration);
p3 = interp1(solubilityParameters(1,:),solubilityParameters(3,:),additiveConcentration);
p2 = interp1(solubilityParameters(1,:),solubilityParameters(4,:),additiveConcentration);
p1 = interp1(solubilityParameters(1,:),solubilityParameters(5,:),additiveConcentration);
p0 = interp1(solubilityParameters(1,:),solubilityParameters(6,:),additiveConcentration);

        % array containing the growth rate parameters for each additive
        % concentration: (c; k1; k2; k3)
growthParameters = [0 0.04 0.08; 0.1 9.7 15; 4e3 5.3e3 5.5e3; 1.48e4 2.04e4 1.8e4];
k1 = interp1(growthParameters(1,:), growthParameters(2,:),additiveConcentration);
k2 = interp1(growthParameters(1,:), growthParameters(3,:),additiveConcentration);
k3 = interp1(growthParameters(1,:), growthParameters(4,:),additiveConcentration);

% equilibrium temperature at initial concentration (starting point of
% any simulation)
T0 = 7.859*fzero(@(theta) 1000*(p4*theta^4 + p3*theta^3 + p2*theta^2 + p1*theta + p0) - initialConcentration,1) + 18.85;

shapeFactor = str2double(essentialParameters(3)); % assuming cuboidal particles
particleDensity = str2double(essentialParameters(4)); % [g/um3]
simulationTime = str2double(essentialParameters(5)); % [h] range of t required

% define length step and spatial domain
dL = str2double(essentialParameters(6)); % [um]
L = 1:dL:1000; % [um]

% prompt user for initial PSD
choicePSD = listdlg('ListString',{'Gaussian','Pulse'},'PromptString','Please select an inital particle size distribution','SelectionMode','single');
% use choice to generate PSD
if choicePSD == 1
    initialPSD = 226458*normpdf(L,300,20);
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
operationMode = listdlg("ListString",{'constant temperature (crash cooling)'},"PromptString",'Please select the required simulation',"SelectionMode","single");
switch operationMode
%% T. Vetter growth rate/solubility (additive)(5) constant temperature comparison
    case 1
        % prompt user for constant temperature and additive concentration
        prompt = {['Please input the constant temperature at which crystallisation takes place [' char(176) 'C]']};
        dlgtitle = 'Constant temperature';
        definput = {'15'};
        answer = inputdlg(prompt,dlgtitle, [1 35],definput);
        constantTemperature = str2double(answer(1));

        % generate temperature ramp (starting from equilibrium)
        temperatureRamp = [0 900 900.01 simulationTime; T0 T0 constantTemperature ...
            constantTemperature];
        
        % high resolution solution:
        [f, concentration, G, supersaturation, m3, t,...
            temperature] = highRes1D_additive(dL, L, simulationTime, k1, k2, k3, p0, p1, p2, p3, p4, shapeFactor,...
            temperatureRamp, particleDensity, initialConcentration, initialPSD, operationMode);
        %calculate average length for comparison to method of moments
        m0 = trapz(L,f);
        m1 = trapz(L,L'.*f);
        averageLength = m1./m0;

        % method of moments solution:
        [t_mom, y] = ode15s(@(t_mom, y)mom_additive(t_mom, y, temperatureRamp, k1, k2, k3, p0, p1, p2, p3, p4, ...
            particleDensity, shapeFactor, operationMode), [0 simulationTime], y0);
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
        
        % plot predicted by method of moments
    
         plot([T0; temperature_mom],[initialConcentration; concentration_mom],'--','linewidth',1.2)
        
        % generate solubility curve for reference
        T = 5:0.1:35;
        theta = (T-18.85)/7.859;
        solubility = 1000*(p4*theta.^4 + p3*theta.^3 + p2*theta.^2 + p1*theta + p0);
        
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
        xlabel({'time' '[s]'})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution','Method of moments')
        
        % average length profile
        subplot(2,1,2)
        % average length profile of system starting at equilibrium then crash
        % cooling:
        plot(t,averageLength,'linewidth',1.2)
        hold on
        
        % average length profile predicted by method of moments
        plot(t_mom,averageLength_mom,'--','linewidth',1.2)
    
        set(gca,'FontSize',18)  
        title('Average length profile')
        xlabel({'time' '[s]'})
        ylabel({'Average length', ['[' char(181) 'm]']})
        legend('High resolution (starting from equilibrium)','Method of moments')
% PSD plots:
        figure(3)
        % PSD evolution comparison
        plot(L,f(:,1),'linewidth',1.2)
        hold on
        plot(L,f(:,end),'linewidth',1.2)
        
        set(gca,'FontSize',18)  
        title('PSD evolution comparison')
        xlabel({'Length' ['[' char(181) 'm]']})
        ylabel({'Number density', '[μm^{-1} kg^{-1}]'})
        legend('Initial PSD (starting from equilibrium)', ...
            'Final PSD (starting from equilibrium)')
end