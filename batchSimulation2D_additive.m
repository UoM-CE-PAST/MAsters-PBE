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
% - 2023/03/20, MA: replaced trapz with sum (more realistic)
% - 2023/03/24, MA: improved high resolution memory efficiency
% - 2024/03/27, MA: added constant supersaturation operation
%
% Purpose: to provide a framework for different simulations of a 2D batch
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
clear;clear mom2D_additiveCS; clc; close all

% create a dialogue box to gather essential parameters from user
prompt = {'Initial concentration [g/kg]','additive concentration [kg/kg]s','shape factor', ['particle density [g/' char(181) 'm' char(179) ']'],'simulation time [h]', ['radius (L1) step size [' char(181) 'm]'], ['height (L2) step size [' char(181) 'm]']};
dlgtitle = 'Essential parameters';
definput = {'8','0','0.5236','1.1128e-12','1000','5','5'};
essentialParameters = inputdlg(prompt,dlgtitle, [1 35],definput);

initialConcentration = str2double(essentialParameters(1)); % [g/kg]

additiveConcentration = str2double(essentialParameters(2)); % [kg/kg]

% additive factor, depends on concentration of additive used
additiveFactors = [1 1.21 1.43;
    1 0.415 0.366];
growthFactor = 1;
solubilityFactor = 1;

kg11 = 3600*58;
kg12 = 2400;
kg13 = 2.5;
kg21 = 3600*2700;
kg22 = 2400;
kg23 = 3.7;

% dissolution rate parameters (independent of additive concentration):
kd11 = 3600*0.272e6;
kd12 = 3223;
kd21 = 3600*1.636e6;
kd22 = 3572;

% equilibrium temperature at initial concentration (starting point of
% any simulation)
T0 = (1/0.036)*log(initialConcentration/(3.37*solubilityFactor));

shapeFactor = str2double(essentialParameters(3)); % assuming cuboidal particles
particleDensity = str2double(essentialParameters(4)); % [g/um3]
simulationTime = str2double(essentialParameters(5)); % [h] range of t required

% define length step and spatial domain
dL1 = str2double(essentialParameters(6)); % [um]
L1 = 1:dL1:800; % [um]
dL2 = str2double(essentialParameters(7));
L2 = 1:dL2:2000; % [um]
[l1, l2] = meshgrid(L1,L2);

L = [l1(:) l2(:)];
clear l1 l2
mean = [200 400];
standardDeviation = [50 0; 0 50];

initialPSD = 1e5*mvnpdf(L,mean, standardDeviation);
initialPSD = reshape(initialPSD,length(L2),length(L1));

% calculate initial moments using intergation of initial PSD
initialm00 = sum(initialPSD,'all')*dL1*dL2;
initialm10 = sum(L1.*initialPSD,'all')*dL1*dL2;
initialm01 = sum(L2'.*initialPSD,'all')*dL1*dL2;
initialm20 = sum(L1.^2.*initialPSD,'all')*dL1*dL2;
initialm30 = sum(L1.^3.*initialPSD,'all')*dL1*dL2;
initialm02 = sum(L2'.^2.*initialPSD,'all')*dL1*dL2;
initialm11 = sum(L1.*L2'.*initialPSD,'all')*dL1*dL2;
initialm21 = sum(L1.^2.*L2'.*initialPSD,'all')*dL1*dL2;
initialm12 = sum(L1.*L2'.^2.*initialPSD,'all')*dL1*dL2;
initialm31 = sum(L1.^3.*L2'.*initialPSD,'all')*dL1*dL2;
initialm22 = sum(L1.^2.*L2'.^2.*initialPSD,'all')*dL1*dL2;



% define vector containing initial conditions for method of moments
y0 = [initialm00 initialm10 initialm01 initialm20 initialm02 initialm11 initialm21 initialm12 initialm22 initialm30 initialm31 initialConcentration];

% prompt user for simulation mode
simulationMode = listdlg("ListString",{'volume-limited temperature cycling','Additive effect on temperature cycles'},"PromptString",'Please select the required simulation',"SelectionMode","single");
switch simulationMode
    %% T. Vetter growth rate/solubility (additive)(5) isovolume operation
    case 1
        % prompt user for desired volume
        prompt = {'Desired final yield (as a multiple of the original mass):','Number of growth/dissolution cycles:','Minimum supersaturation:','Maximum supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'5','2','0.9','1.15'};
        cycleParameters = inputdlg(prompt,dlgtitle, [1 35],definput);
        yieldFactor = str2double(cycleParameters(1));
        maxCycleNumber = str2double(cycleParameters(2));
        Smin = str2double(cycleParameters(3));
        Smax = str2double(cycleParameters(4));
        supersaturationLimits = [Smin Smax];
        m21min = initialm21;
        m21max = yieldFactor*m21min;
        cmax = initialConcentration;
        cmin = cmax - particleDensity*shapeFactor*(m21max - m21min);
        Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
        Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

        % memory efficient high resolution solution
        [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
        averageRadius = m31./m21;
        averageHeight = m22./m21;
        % method of moments solution:
        [t_mom, y] = ode15s(@(t_mom, y)mom2D_additiveCS(t_mom,y,kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,shapeFactor,particleDensity,initialConcentration,initialm21,yieldFactor,maxCycleNumber, supersaturationLimits, growthFactor,solubilityFactor), [0 simulationTime], y0);
        concentration_mom = y(:,12);
        m00_mom = y(:,1);
        m21_mom = y(:,7);
        m31_mom = y(:,11);
        m22_mom = y(:,9);
        averageRadius_mom = m31_mom./m21_mom;
        averageHeight_mom = m22_mom./m21_mom;
        % extract temperature profile used by method of moments:
        temperature_mom = zeros(length(y),1);
        for i = 1:length(y)
            [~, temperature_mom(i)] = mom2D_additiveCS(t_mom(i),y(i,:),kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,shapeFactor,particleDensity,initialConcentration,initialm21,yieldFactor,maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
        end

        % aspect ratio plots:
        figure(1)
        plot(averageRadius, averageHeight,'LineWidth',1.2)
        hold on
        plot(averageRadius_mom, averageHeight_mom, '--','LineWidth',1.2)
        % ideal aspect ratio:
        plot([0 xlim],[0 xlim],':','LineWidth',1.2)
        % iso-volume plots:
        l1 = 200:0.1:550;
        avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
        l2 = avgVolumemin./l1.^2;
        plot(l1,l2,'LineWidth',1.2) % lower bound
        l2 = yieldFactor*avgVolumemin./l1.^2;
        plot(l1,l2,'LineWidth',1.2) % desired final volume
        set(gca,'FontSize',18)
        title('Aspect ratio')
        xlabel({'L1', ['[' char(181) 'm]']})
        ylabel({'L2', ['[' char(181) 'm]']})
        legend('High resolution','Method of moments','Aspect ratio = 1','Minimum volume','Desired volume')

        % concentration plot
        figure(2)
        plot(t,concentration,'linewidth',1.2)
        hold on
        plot(t_mom,concentration_mom,'--','linewidth',1.2)
        set(gca,'FontSize',18)
        title('Concentration profile')
        xlabel({'time' '[s]'})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution','Method of moments')

        % solubility plot
        figure(3)
        plot(temperature,concentration,'linewidth',1.6)
        hold on
        plot([T0; temperature_mom],[initialConcentration; concentration_mom],'--','linewidth',1.2)
        % generate solubility curve for reference
        T = 0:0.1:40;
        solubility =solubilityFactor*3.37*exp(0.036*T);
        plot(T,solubility,'linewidth',1.2)
        set(gca,'FontSize',18)
        title('Solubility curve')
        xlabel({'Temperature', ['[' char(176) 'C]']})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution', ...
            'Method of moments ', ...
            'Solubility curve')
        xline([Tmin Tmax],'--',{'Minimum temperature','Maximum temperature'})

        % PSSD plots
        figure(4)
        contourf(L1,L2,initialPSD,[5 100 200 300 400]) % initial PSSD
        hold on
        contourf(L1,L2,finalPSD, [5 100 200 300 400],'--') % final PSSD
        xlim([100 600])
        ylim([200 1400])

        set(gca,'FontSize',18)
        title('PSSD')
        xlabel({'L1', ['[' char(181) 'm]']})
        ylabel({'L2', ['[' char(181) 'm]']})

        %         figure(5)
        %         surf(L1,L2,initialPSD) % initial PSSD
        %         hold on
        %         surf(L1,L2,finalPSD) % final PSSD
        %         xlim([100 400])
        %         ylim([200 800])
        %
        %         set(gca,'FontSize',18)
        %         title('PSSD')
        %         xlabel({'L1', ['[' char(181) 'm]']})
        %         ylabel({'L2', ['[' char(181) 'm]']})

        %% Additive concentration effect on temperature cycles
    case 2
        prompt = {'Desired final yield (as a multiple of the original mass):','Number of growth/dissolution cycles:','Minimum supersaturation:','Maximum supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'20','2','0.9','1.15'};
        cycleParameters = inputdlg(prompt,dlgtitle, [1 35],definput);
        yieldFactor = str2double(cycleParameters(1));
        maxCycleNumber = str2double(cycleParameters(2));
        Smin = str2double(cycleParameters(3));
        Smax = str2double(cycleParameters(4));
        supersaturationLimits = [Smin Smax];
        m21min = initialm21;
        m21max = yieldFactor*initialm21;
        cmax = initialConcentration;
        cmin = cmax - particleDensity*shapeFactor*(m21max - m21min);


        for k = 1:length(additiveFactors)
            growthFactor = additiveFactors(2,k); % slows down rate of radius growth appropriately
            p4 = solubilityParameters(2,k);
            p3 = solubilityParameters(3,k);
            p2 = solubilityParameters(4,k);
            p1 = solubilityParameters(5,k);
            p0 = solubilityParameters(6,k);
            Tmax = 7.859*fzero(@(theta) 1000*(p4*theta^4 + p3*theta^3 ...
                + p2*theta^2 + p1*theta + p0) - cmax,1) + 18.85;
            Tmin = 7.859*fzero(@(theta) 1000*(p4*theta^4 + p3*theta^3 ...
                + p2*theta^2 + p1*theta + p0) - cmin,1) + 18.85;
            % memory efficient high resolution solution
        [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, p0, p1, p2, p3, p4, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor);
        averageRadius = m31./m21;
        averageHeight = m22./m21;
        % method of moments solution:
        [t_mom, y] = ode15s(@(t_mom, y)mom2D_additiveCS(t_mom,y,kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,p0,p1,p2,p3,p4,shapeFactor,particleDensity,initialConcentration,initialm21,yieldFactor,maxCycleNumber, supersaturationLimits, growthFactor), [0 simulationTime], y0);
        concentration_mom = y(:,12);
        m00_mom = y(:,1);
        m21_mom = y(:,7);
        m31_mom = y(:,11);
        m22_mom = y(:,9);
        averageRadius_mom = m31_mom./m21_mom;
        averageHeight_mom = m22_mom./m21_mom;
        % extract temperature profile used by method of moments:
        temperature_mom = zeros(length(y),1);
        for i = 1:length(y)
            [~, temperature_mom(i)] = mom2D_additiveCS(t_mom(i),y(i,:),kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,p0,p1,p2,p3,p4,shapeFactor,particleDensity,initialConcentration,initialm21,yieldFactor,maxCycleNumber, supersaturationLimits,growthFactor);
        end

        % aspect ratio plots:
        figure(1)
        plot(averageRadius, averageHeight,'LineWidth',1.2)
        hold on
        plot(averageRadius_mom, averageHeight_mom, '--','LineWidth',1.2)
        % ideal aspect ratio:
        plot([0 xlim(2)],[0 xlim(2)],':','LineWidth',1.2)
        % iso-volume plots:
        l1 = 200:0.1:550;
        avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
        l2 = avgVolumemin./l1.^2;
        plot(l1,l2,'LineWidth',1.2) % lower bound
        l2 = yieldFactor*avgVolumemin./l1.^2;
        plot(l1,l2,'LineWidth',1.2) % desired final volume
        set(gca,'FontSize',18)
        title('Aspect ratio')
        xlabel({'L1', ['[' char(181) 'm]']})
        ylabel({'L2', ['[' char(181) 'm]']})
        legend('High resolution','Method of moments','Aspect ratio = 1','Minimum volume','Desired volume')

        % concentration plot
        figure(2)
        plot(t,concentration,'linewidth',1.2)
        hold on
        plot(t_mom,concentration_mom,'--','linewidth',1.2)
        set(gca,'FontSize',18)
        title('Concentration profile')
        xlabel({'time' '[s]'})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution','Method of moments')

        % solubility plot
        figure(3)
        plot(temperature,concentration,'linewidth',1.6)
        hold on
        plot([T0; temperature_mom],[initialConcentration; concentration_mom],'--','linewidth',1.2)
        % generate solubility curve for reference
        T = 5:0.1:35;
        solubility = solubilityFactor*3.37*exp(0.036*T);
        plot(T,solubility,'linewidth',1.2)
        set(gca,'FontSize',18)
        title('Solubility curve')
        xlabel({'Temperature', ['[' char(176) 'C]']})
        ylabel({'Concentration' '[g kg^{-1}]'})
        legend('High resolution', ...
            'Method of moments ', ...
            'Solubility curve')
        xline([Tmin Tmax],'--',{'Minimum temperature','Maximum temperature'})

        % PSSD plots
        figure(4)
        contourf(L1,L2,initialPSD,[5 100 200 300 400]) % initial PSSD
        hold on
        contourf(L1,L2,finalPSD, [5 100 200 300 400],'--') % final PSSD
        xlim([100 600])
        ylim([200 1400])

        set(gca,'FontSize',18)
        title('PSSD')
        xlabel({'L1', ['[' char(181) 'm]']})
        ylabel({'L2', ['[' char(181) 'm]']})
        end
end