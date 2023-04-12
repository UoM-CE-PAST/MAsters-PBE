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
% - 2023/04/01, MA: added 1D numerical scheme comparison
% - 2023/04/02, MA: improved plots + added automatic plot saving
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
prompt = {'Initial concentration [g/kg]','additive concentration [kg/kg]','shape factor', ['particle density [g/' char(181) 'm' char(179) ']'],'simulation time [h]', ['length step size [' char(181) 'm]']};
dlgtitle = 'Essential parameters';
definput = {'8','0','0.5236','1.1128e-12','50','1'};
essentialParameters = inputdlg(prompt,dlgtitle, [1 35], definput);

initialConcentration = str2double(essentialParameters(1)); % [g/kg]

additiveConcentration = str2double(essentialParameters(2)); % [kg/kg]

% additive effect (NEEDS UPDATE)
solubilityFactor = 1;
growthFactor = 1;

% kinetics
kg1 = 3600*2700;
kg2 = 2400;
kg3 = 3.7;
kd1 = 3600*1.636e6/1000;
kd2 = 3572;

% equilibrium temperature at initial concentration (starting point of
% any simulation)
T0 = (1/0.036)*log(initialConcentration/(3.37*solubilityFactor));

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
    initialPSD = 1e5*normpdf(L,300,20);
elseif choicePSD == 2
    initialPSD = zeros(1,length(L));
    initialPSD(250:350)=1e3;
end

% calculate initial moments using intergation
m0i=trapz(L,initialPSD);
m1i=trapz(L,L.*initialPSD);
m2i=trapz(L,L.^2.*initialPSD);
m3i=trapz(L,L.^3.*initialPSD);
m4i = trapz(L,L.^4.*initialPSD);
m5i = trapz(L,L.^5.*initialPSD);

% define vector containing initial conditions for method of moments
y0 = [m0i m1i m2i m3i m4i m5i initialConcentration];

% prompt user for operation mode
simulationMode = listdlg("ListString",{'constant temperature (crash cooling)','Numerical scheme comparison'},"PromptString",'Please select the required simulation',"SelectionMode","single");
switch simulationMode
    %% T. Vetter growth rate/solubility (additive)(5) constant temperature comparison
    case 1
        % prompt user for constant temperature
        prompt = {['Please input the constant temperature at which crystallisation takes place [' char(176) 'C]']};
        dlgtitle = 'Constant temperature';
        definput = {'10'};
        answer = inputdlg(prompt,dlgtitle, [1 35],definput);
        constantTemperature = str2double(answer(1));

        % generate temperature ramp (starting from equilibrium)
        temperatureRamp = [0 1 1.001 simulationTime; T0 T0 constantTemperature ...
            constantTemperature];

        % high resolution solution:
        [f, concentration, G, supersaturation, m3, t,...
            temperature] = highRes1D_additiveTD(dL, L, simulationTime, kg1, kg2, kg3, kd1, kd2, shapeFactor,...
            temperatureRamp, particleDensity, initialConcentration, initialPSD, growthFactor, solubilityFactor);
        %calculate average length for comparison to method of moments
        m0 = trapz(L,f);
        m1 = trapz(L,L'.*f);
        m4 = trapz(L,L'.^4.*f);
        m5 = trapz(L,L'.^5.*f);
        averageLength = m4./m3;
        stdLength = sqrt(abs(m5./m3-averageLength.^2));

        % method of moments solution:
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t_mom, y] = ode113(@(t_mom, y)mom_additive(t_mom, y, temperatureRamp, kg1, kg2, kg3, kd1, kd2, particleDensity, ...
            shapeFactor, growthFactor, solubilityFactor), [0 simulationTime], y0,options);
        concentration_mom = y(:,7);
        % extract temperature profile predicted by method of moments
        temperature_mom = interp1(temperatureRamp(1,:),temperatureRamp(2,:),t_mom);
        % average length predicted by method of moments
        averageLength_mom  = y(:,5)./y(:,4);
        stdLength_mom = sqrt(abs(y(:,6)./y(:,4)-averageLength_mom.^2));

        % solubility curve plots:
        figure(1)
        % plot of system starting at equilibrium then crash cooling
        plot(temperature,concentration,'linewidth',1.6)
        hold on

        % plot predicted by method of moments

        plot(temperature_mom,concentration_mom,'--','linewidth',1.2)

        % generate solubility curve for reference
        T = 0:0.1:45;

        solubility = solubilityFactor*3.37*exp(0.036*T);

        % plot of solubility curve to show driving force
        plot(T,solubility,'linewidth',1.2)

        set(gca,'FontSize',18)
        title('Solubility curve')
        xlabel({'Temperature', ['[' char(176) 'C]']})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution', ...
            'Method of moments ', ...
            'Solubility curve')

        % concentration and average length plots:
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

        % average length profile predicted by method of moments
        plot(t_mom,averageLength_mom,'--','linewidth',1.2)

        set(gca,'FontSize',18)
        title('Average length profile')
        xlabel({'time' '[h]'})
        ylabel({'Average length', ['[' char(181) 'm]']})
        legend('High resolution (starting from equilibrium)','Method of moments')
        % PSD plots:
        figure(3)
        % PSD evolution comparison
        plot(L,f(:,1),'linewidth',1.2)
        hold on
        plot(L,f(:,end),'linewidth',1.2)

        set(gca,'FontSize',18)
        title('PSD evolution')
        xlabel({'Length' ['[' char(181) 'm]']})
        ylabel({'Number density', '[μm^{-1} kg^{-1}]'})
        legend('Initial PSD (starting from equilibrium)', ...
            'Final PSD (starting from equilibrium)')
        % standard deviation plots:
        figure(4)
        plot(t,stdLength,'LineWidth',1.2)
        hold on
        plot(t_mom,stdLength_mom,'--','LineWidth',1.2)

        set(gca,'FontSize',18)
        title('Standard deviation profile')
        xlabel({'time', '[h]'})
        ylabel({'Number density', '[μm^{-1} kg^{-1}]'})
        legend('High resolution', ...
            'Method of moments')
        % moments plots
        figure(5)
        subplot(3,2,1)
        plot(t,m0,'LineWidth',1.2)
        hold on
        plot(t_mom,y(:,1),'--','LineWidth',1.2)
        set(gca,'FontSize',18)
        title('m0 profile')

        subplot(3,2,2)
        plot(t,m1,'LineWidth',1.2)
        hold on
        plot(t_mom,y(:,2),'--','LineWidth',1.2)
        set(gca,'FontSize',18)
        title('m1 profile')

        subplot(3,2,3)
        plot(t,m4,'LineWidth',1.2)
        hold on
        plot(t_mom,y(:,5),'--','LineWidth',1.2)
        set(gca,'FontSize',18)
        title('m4 profile')

        subplot(3,2,4)
        plot(t,m5,'LineWidth',1.2)
        hold on
        plot(t_mom,y(:,6),'--','LineWidth',1.2)
        set(gca,'FontSize',18)
        title('m5 profile')

        subplot(3,2,[5 6])
        plot(t,m3,'LineWidth',1.2)
        hold on
        plot(t_mom,y(:,4),'--','LineWidth',1.2)
        set(gca,'FontSize',18)
        title('m3 profile')

        %% Numerical schemes comparison
    case 2
% create a folder with the current date and time:
        timeDate = string(datetime('now','TimeZone','local','Format', 'dd-MM-yyyy HH.mm.ss'));
        mkdir('C:\Users\moham\OneDrive - The University of Manchester\Year 4 modules\Dissertation (CHEN40100)\Figures\Batch1D_additive\SchemeComparison', timeDate)
% prompt user for constant temperature
        prompt = {['Please input the constant temperature at which crystallisation takes place [' char(176) 'C]']};
        dlgtitle = 'Constant temperature';
        definput = {'-10'};
        answer = inputdlg(prompt,dlgtitle, [1 35],definput);
        constantTemperature = str2double(answer(1));

        % generate temperature ramp (starting from equilibrium)
        temperatureRamp = [0 simulationTime; constantTemperature ...
            constantTemperature];
% upwind comparison:
        % Gaussian:
        initialPSD = 1e5*normpdf(L,300,20);

% high resolution solution:
        [f, ~, ~, ~, m3, t] = highRes1D_additiveTD(dL, L, simulationTime, kg1, kg2, kg3, kd1, kd2, shapeFactor,...
            temperatureRamp, particleDensity, initialConcentration, initialPSD, growthFactor, solubilityFactor);
        %calculate average length for comparison to method of moments
        m0 = trapz(L,f);
        m1 = trapz(L,L'.*f);
        m4 = trapz(L,L'.^4.*f);
        m5 = trapz(L,L'.^5.*f);
        averageLength = m4./m3;
        stdLength = sqrt(abs(m5./m3-averageLength.^2));

% upwind solution:
        [f_upwind, ~, ~, ~, m3_upwind, t_upwind] = upwind_additiveTD(dL, L, simulationTime, kg1, kg2, kg3, kd1, kd2, shapeFactor,...
            temperatureRamp, particleDensity, initialConcentration, initialPSD, growthFactor, solubilityFactor);
        %calculate average length for comparison to method of moments
        m0_upwind = trapz(L,f);
        m1_upwind = trapz(L,L'.*f);
        m4_upwind = trapz(L,L'.^4.*f);
        m5_upwind = trapz(L,L'.^5.*f);
        averageLength_upwind = m4./m3_upwind;
        stdLength_upwind = sqrt(abs(m5./m3_upwind-averageLength_upwind.^2));

% standard deviation comparison:
% Lax-Wendroff solution (other standard deviations already available):
        [~, ~, ~, ~, m3_LW, t_LW] = LaxWendroff_additiveTD(dL, L, simulationTime, kg1, kg2, kg3, kd1, kd2, shapeFactor,...
            temperatureRamp, particleDensity, initialConcentration, initialPSD, growthFactor, solubilityFactor);
        m0_LW = trapz(L,f);
        m1_LW = trapz(L,L'.*f);
        m4_LW = trapz(L,L'.^4.*f);
        m5_LW = trapz(L,L'.^5.*f);
        averageLength_LW = m4./m3_LW;
        stdLength_LW = sqrt(abs(m5./m3_LW-averageLength_LW.^2));

% standard deviation plots:
        figure(1)
        % fraction of simulation time to stop plotting:
        upperlimit = simulationTime; 
        lowerlimit = 0.1*simulationTime;

        plot(t,stdLength/stdLength(1),'LineWidth',1.5)
        hold on

        plot(t_upwind,stdLength_upwind/stdLength(1),'-.','LineWidth',1.5)

        plot(t_LW,stdLength_LW/stdLength(1),'--','LineWidth',1.5)

        set(gca,'FontSize',8,'XTick',[],'fontname','times')
        fig = gcf;
        fig.Position(3:4)=[312 312];
        colororder(['#000000';'#00FFFF';'#FFD700'])
        legend('High resolution scheme', 'Upwind scheme', 'Lax-Wendroff scheme','Location','northwest')
        xlabel({'{\it t}','[h]'})
        ylabel({'\sigma_{model}/\sigma_{actual}','[-]'})
        xlim([lowerlimit upperlimit])
        fileName= 'Gaussian standard deviation scheme comaprison.pdf';
        fpath = strcat('C:\Users\moham\OneDrive - The University of Manchester\Year 4 modules\Dissertation (CHEN40100)\Figures\Batch1D_additive\SchemeComparison\', timeDate);
        exportgraphics(gca, fullfile(fpath, fileName), 'ContentType','vector');
        clear upperlimit lower limit

% upwind PSD comparison plots:
        figure(2)
        plot(L, initialPSD,':','LineWidth',1.5) % initial
        hold on

        plot(L,f(:,end), 'LineWidth',1.5) % high resolution

        plot(L, f_upwind(:,end),'--', 'LineWidth',1.5) % upwind

        set(gca,'FontSize',8,'FontName','times')
        fig = gcf;
        fig.Position(3:4)=[312 312];
        colororder(["#000000";"#000000";"#FFD700"])
        legend('Initial condition', 'High resolution solution', 'Upwind solution','Location','north')
        xlabel({'L','[\mum]'})
        ylabel({'f','[kg^{-1}\mum^{-1}]'})

% Lax-Wendroff comparison:
        % pulse:
        initialPSD = zeros(1,length(L));
        initialPSD(250:350)=1e3; 

% high resolution solution:
        f = highRes1D_additiveTD(dL, L, simulationTime, kg1, kg2, kg3, kd1, kd2, shapeFactor,...
            temperatureRamp, particleDensity, initialConcentration, initialPSD, growthFactor, solubilityFactor);

% Lax-Wendroff solution:
        f_LW = LaxWendroff_additiveTD(dL, L, simulationTime, kg1, kg2, kg3, kd1, kd2, shapeFactor,...
            temperatureRamp, particleDensity, initialConcentration, initialPSD, growthFactor, solubilityFactor);
        
% Lax-Wendroff PSD comparison plots:
        figure(3)
        plot(L, initialPSD,':','LineWidth',1.5) % initial
        hold on

        plot(L,f(:,end), 'LineWidth',1.5) % high resolution

        plot(L, f_LW(:,end),'--', 'LineWidth',1.5) % Lax-Wendroff

        set(gca,'FontSize',8,'fontname','times')
        fig = gcf;
        fig.Position(3:4)=[312 312];
        colororder(['#000000';'#000000';'#00A300'])
        legend('Initial condition', 'High resolution solution', 'Lax-Wendroff solution','Location','north')
        xlabel({'{\it L}',['[' char(181) 'm]']})
        ylabel({'{\it f}',['[kg^{-1}' char(181) 'm^{-1}]']})
        fileName='Lax-Wendroff PSD comparison.pdf';
        fpath = strcat('C:\Users\moham\OneDrive - The University of Manchester\Year 4 modules\Dissertation (CHEN40100)\Figures\Batch1D_additive\SchemeComparison\', timeDate);
        exportgraphics(gca, fullfile(fpath, fileName), 'ContentType','vector');

% store operating conditions etc in a matrix to be exported as a
% text file
        namesCell = {'Initial Concentration'; 'Additive concentration';'Shape factor'; 'Particle density'; 'Simulation time';'Bin size';'Temperature'};
        unitsCell = {'[g/kg]';'[kg/kg]'; '[-]';['[g/' char(181) 'm' char(179) ']'];'[h]';['[' char(181) 'm]'];['[' char(176) 'C]']};
        parametersCell = cat(1, [namesCell, [essentialParameters(:); constantTemperature], unitsCell]);
        fileName='Numerical schemes comparison parameters';
        fpath = strcat('C:\Users\moham\OneDrive - The University of Manchester\Year 4 modules\Dissertation (CHEN40100)\Figures\Batch1D_additive\SchemeComparison\', timeDate);
        writecell(parametersCell,fullfile(fpath,fileName))
end