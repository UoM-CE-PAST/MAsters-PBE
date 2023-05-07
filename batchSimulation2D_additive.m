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
% - 2023/03/27, MA: added constant supersaturation operation
% - 2023/04/23, MA: massive improvements to simulation and plots
% - 2023/05/06, MA: added attainable region analysis
% - 2023/05/07, MA: major clean up
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
clear; clc; close all
[status,hash] = system(strcat(['git rev-list -1 HEAD ' mfilename '.m']));
hash = hash(1:7);
if status ~= 0
    hash = [];
    warning('git does not recognise current file. Check if the file has been commited.')
end
timeDate = string(datetime('now','TimeZone','local','Format', 'dd-MM-yyyy HH.mm.ss'));
folderName = strcat(timeDate,{' '},hash);

% path for figures to be saved in, change if required:
path = 'C:\Users\moham\OneDrive - The University of Manchester\Year 4 modules\Dissertation (CHEN40100)\Figures\';

% create a dialogue box to gather essential parameters from user:
prompt = {'Initial concentration [g/kg]','additive concentration [kg/kg]s','shape factor', ['particle density [g/' char(181) 'm' char(179) ']'],'simulation time [h]', ['radius (L1) step size [' char(181) 'm]'], ['height (L2) step size [' char(181) 'm]']};
dlgtitle = 'Essential parameters';
definput = {'8','0','0.7854','1.1128e-12','1000','1','1'};
essentialParameters = inputdlg(prompt,dlgtitle, [1 35],definput);

initialConcentration = str2double(essentialParameters(1)); % [g/kg]

additiveConcentration = str2double(essentialParameters(2)); % [kg/kg]

% additive factor, depends on concentration of additive used
additiveFactors = [0 0.04 0.08;
    1 1.21 1.43; % solubility factors
    1 0.415 0.366]; % growth rate factors

solubilityFactor = interp1(additiveFactors(1,:),additiveFactors(2,:),additiveConcentration);
growthFactor = interp1(additiveFactors(1,:),additiveFactors(3,:),additiveConcentration);

kg11 = 3600*58;
kg12 = 2400;
kg13 = 2.5;
kg21 = 3600*0.5*2700;
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
L1 = 1:dL1:600; % [um]
dL2 = str2double(essentialParameters(7));
L2 = 1:dL2:1000; % [um]
[l1, l2] = meshgrid(L1,L2);

L = [l1(:) l2(:)];
clear l1 l2
meanDimensions = [250 400];
standardDeviation = [50 0; 0 50];

initialPSD = 1e5*mvnpdf(L,meanDimensions, standardDeviation);
initialPSD = reshape(initialPSD,length(L2),length(L1));

% calculate initial moments using intergation of initial PSD
initialm00 = sum(initialPSD,'all')*dL1*dL2;

initialm10 = sum(L1.*initialPSD,'all')*dL1*dL2;
initialm01 = sum(L2'.*initialPSD,'all')*dL1*dL2;
initialm40 = sum(L1.^4.*initialPSD,'all')*dL1*dL2;
initialm30 = sum(L1.^3.*initialPSD,'all')*dL1*dL2;
initialm20 = sum(L1.^2.*initialPSD,'all')*dL1*dL2;
initialm02 = sum(L2'.^2.*initialPSD,'all')*dL1*dL2;
initialm03 = sum(L2'.^3.*initialPSD,'all')*dL1*dL2;
initialm11 = sum(L1.*L2'.*initialPSD,'all')*dL1*dL2;
initialm21 = sum(L1.^2.*L2'.*initialPSD,'all')*dL1*dL2;
initialm12 = sum(L1.*L2'.^2.*initialPSD,'all')*dL1*dL2;
initialm31 = sum(L1.^3.*L2'.*initialPSD,'all')*dL1*dL2;
initialm13 = sum(L1.*L2'.^3.*initialPSD,'all')*dL1*dL2;
initialm41 = sum(L1.^4.*L2'.*initialPSD,'all')*dL1*dL2;
initialm22 = sum(L1.^2.*L2'.^2.*initialPSD,'all')*dL1*dL2;
initialm23 = sum(L1.^2.*L2'.^3.*initialPSD,'all')*dL1*dL2;

% define vector containing initial conditions for method of moments
y0 = [initialm00 initialm10 initialm01 initialm20 initialm02 initialm11 initialm21 initialm12 initialm22 initialm30 initialm31 initialm03 initialm40 initialm13 initialm23 initialm41 initialConcentration];

% prompt user for simulation mode
simulationMode = listdlg("ListString",{'volume-limited temperature cycling','Additive effect on temperature cycles','Supersaturation limit effect on temperature cycles','Relative kinetics effect on temperature cycles','Pure growth base case','Pure growth bin size varitaion','Attainable region pure growth','Powerpoint figures'},"PromptString",'Please select the required simulation',"SelectionMode","single");
switch simulationMode
    %% Iso-volume temperature cycles base case
    case 1
% create a folder with the current date and time:
        mkdir(strcat(path,'Batch2D_additive\Iso-volume base case cycles'), folderName)

% prompt user for desired volume and supersaturation limits
        prompt = {'Desired final yield (as a multiple of the original mass):','Number of growth/dissolution cycles:','Minimum supersaturation:','Maximum supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'3','2','0.9','1.15'};
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

% memory efficient high resolution solution:
        [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
        averageRadius = m31./m21;
        averageHeight = m22./m21;
        stdRadius = sqrt(abs(m41./m21 - averageRadius.^2));
        stdHeight = sqrt(abs(m23./m21 - averageHeight.^2));
        temperatureRamp = [t; temperature]; % to be used by mom
        
% method of moments solution:
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t_mom, y] = ode113(@(t_mom, y)mom2D_additiveTD(t_mom,y,kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,temperatureRamp,shapeFactor,particleDensity,solubilityFactor,growthFactor), [0 simulationTime], y0,options);
        concentration_mom = y(:,17);
        m00_mom = y(:,1);
        m21_mom = y(:,7);
        m31_mom = y(:,11);
        m22_mom = y(:,9);
        m41_mom = y(:,16);
        m23_mom = y(:,15);
        averageRadius_mom = m31_mom./m21_mom;
        averageHeight_mom = m22_mom./m21_mom;
        stdRadius_mom = sqrt(abs(m41_mom./m21_mom - averageRadius_mom.^2));
        stdHeight_mom = sqrt(abs(m23_mom./m21_mom - averageHeight_mom.^2));

% aspect ratio plots:
        figure(1)
        tiles1 = tiledlayout(1,2,'Padding','compact','TileSpacing','compact','Units','inches');
        tiles1.InnerPosition(1:2) = [3 3];
        tiles1.InnerPosition(4) = 3.25;
        tiles1.OuterPosition(3) = 7;
        nexttile
        plot(averageRadius, averageHeight,"LineWidth",1)
        hold on
        plot(averageRadius_mom, averageHeight_mom, '--',"LineWidth",1)
        
        % iso-volume plots:
        l1 = 100:0.1:1000;
        avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
        l2 = avgVolumemin./l1.^2;
        plot(l1,l2,"LineWidth",1) % lower bound
        l2 = yieldFactor*avgVolumemin./l1.^2;
        plot(l1,l2,"LineWidth",1) % desired final volume

        % ideal aspect ratio:
        plot(L2,L2,':',"LineWidth",1)
        xlim([200 600]);
        ylim([200 800]);

        % PSSD plots
        contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.6); % initial PSSD
        contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.6); % final PSSD


        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        colororder(gca,['#000000'; '#00A300';'#d00000';'#03045e'; '#000000'])
        xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',10)
        xticks(linspace(200,600,3))
        ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',10)
        yticks(linspace(200,800,4))
        % pbaspect([1 1.1 1])
        legend('High resolution solution','Method of moments solution','Minimum volume','Desired volume','Ideal aspect ratio','Location','southeast')
        title('a)','FontWeight','bold','Fontsize',12)
        text([averageRadius(1) averageRadius(end)], [averageHeight(1) averageHeight(end)],{'     start' '     end'})
        
% solubility plot:
        figure(1)
        nexttile
        plot(temperature,concentration,"LineWidth",1)
        hold on
        plot(temperature(supersaturation>=1),concentration(supersaturation>=1),"LineWidth",1)
        
        % generate solubility curve for reference
        T = 0:0.1:40;
        solubility =solubilityFactor*3.37*exp(0.036*T);
        plot(T,solubility,"LineWidth",1)
        
        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        xlabel({['{\it T} [' char(176) 'C]']},'FontSize',10)
        ylabel({'{\it c} [g kg^{-1}]'},'FontSize',10)
        % pbaspect([1 1.1 1])
        legend('Dissolution stage','Growth stage', ...
            'Solubility curve','AutoUpdate','off','Location','southeast');
        xline([Tmin Tmax],'--',{'{\it T}_{min}','{\it T}_{max}'})
        title('b)','FontWeight','bold','Fontsize',12)
        colororder(gca,['#D00000';'#023E8A';'#000000'])
        text([temperature(1) temperature(end)], [concentration(1) concentration(end)],{'\leftarrow start' '\leftarrow end'})

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('aspect ratio + solubility base case.pdf');
        fpath = strcat(path,'Batch2D_additive\Iso-volume base case cycles\', folderName);
        exportgraphics(tiles1, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

% concentration plot:
        figure(2)
        tiles2 = tiledlayout(3,1,'Padding','compact','TileSpacing','compact','Units','inches');
        tiles2.InnerPosition(1:2) = [3 3];
        tiles2.InnerPosition(4) = 6;
        tiles2.OuterPosition(3) = 7;
        nexttile
        plot(t,concentration,"LineWidth",1)
        hold on
        plot(t_mom,concentration_mom,'--',"LineWidth",1)

        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        xlabel({'{\it t} [h]'},'FontSize',10)
        ylabel({'{\it c} [g kg^{-1}]'},'FontSize',10)
        title('a)','FontWeight','bold','Fontsize',12)
        colororder(gca,['#000000';'#00A300'])
        leg = legend('High resolution solution', 'Method of moments solution','AutoUpdate','off','Orientation','horizontal','Box','off','FontSize',10);
% standard deviation plots
        nexttile
        plot(t,stdRadius./stdRadius(1),"LineWidth",1)
        hold on
        plot(t_mom,stdRadius_mom./stdRadius_mom(1),'--',"LineWidth",1)

        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        xlabel({'{\it t} [h]'},'FontSize',10)
        ylabel({'\sigma_{11,model}/\sigma_{11,actual} [-]'},'FontSize',10)
        title('b)','FontWeight','bold','Fontsize',12)
        colororder(gca,['#000000';'#00A300'])

        nexttile
        plot(t,stdHeight./stdHeight(1),"LineWidth",1)
        hold on
        plot(t_mom,stdHeight_mom./stdHeight_mom(1),'--',"LineWidth",1)

        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        xlabel({'{\it t} [h]'},'FontSize',10)
        ylabel({'\sigma_{22,model}/\sigma_{22,actual} [-]'},'FontSize',10)
        title('c)','FontWeight','bold','Fontsize',12)
        colororder(gca,['#000000';'#00A300'])
        
        leg.Layout.Tile = 'south';

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];
        
        fileName = strcat('Concentration + standard deviation profiles base case.pdf');
        fpath = strcat(path,'Batch2D_additive\Iso-volume base case cycles\', folderName);
        exportgraphics(tiles2, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

% PSSD plots
        figure(3)
        contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9]) % initial PSSD
        hold on
        contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'--') % final PSSD
        xlim([150 850]);
        ylim([150 850]);

        set(gca,'FontSize',8,'FontName','times')
        fig = gcf;
        fig.Units = "inches";
        fig.OuterPosition(3)=3.25;
        title('PSSD')
        xlabel({'L1', ['[' char(181) 'm]']})
        ylabel({'L2', ['[' char(181) 'm]']})

% moment plots:
        figure(4)
        subplot(3,2,1)
        plot(t,m00,"LineWidth",1)
        hold on
        plot(t_mom,m00_mom,'--',"LineWidth",1)

        set(gca,'FontSize',20)
        title('m00 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,2)
        plot(t,m21,"LineWidth",1)
        hold on
        plot(t_mom,m21_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m21 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,3)
        plot(t,m31,"LineWidth",1)
        hold on
        plot(t_mom,m31_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m31 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,4)
        plot(t,m22,"LineWidth",1)
        hold on
        plot(t_mom,m22_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m22 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,5)
        plot(t,m41,"LineWidth",1)
        hold on
        plot(t_mom,m41_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m41 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,6)
        plot(t,m23,"LineWidth",1)
        hold on
        plot(t_mom,m23_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m23 profile')
        legend('High resolution','Method of moments')

% store operating conditions etc in a matrix to be exported as a
% text file
        namesCell = {'Initial Concentration'; 'Additive concentration';'Shape factor'; 'Particle density'; 'Simulation time';'dL1';'dL2';'Yield';'Number of cycles';'Minimum supersaturation';'Maximum supersaturation'};
        unitsCell = {'[g/kg]';'[kg/kg]'; '[-]';['[g/' char(181) 'm' char(179) ']'];'[h]';['[' char(181) 'm]'];['[' char(181) 'm]'];'[-]';'[-]';'[-]';'[-]'};
        parametersCell = cat(1, [namesCell, [essentialParameters(:); cycleParameters], unitsCell]);
        fileName='Iso-volume base case parameters';
        fpath = strcat(path,'Batch2D_additive\Iso-volume base case cycles\', folderName);
        writecell(parametersCell,fullfile(fpath,fileName))
 
%% Additive concentration effect on temperature cycles
    case 2
% create a folder with the current date and time:
        mkdir(strcat(path,'Batch2D_additive\Additive variation cycles'), folderName)

% prompt user for desired volume and supersaturation limits
        prompt = {'Desired final yield (as a multiple of the original mass):','Number of growth/dissolution cycles:','Minimum supersaturation:','Maximum supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'3','2','0.9','1.15'};
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
        Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
        Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

        additiveRange = [0 0.04 0.08];
        subplotTitles = {'a)', 'b)', 'c)'};
        for k = 1:length(additiveRange)
            additiveConcentration = additiveRange(k);
            solubilityFactor = interp1(additiveFactors(1,:),additiveFactors(2,:),additiveRange(k)); % increase solubility appropriately
            growthFactor = interp1(additiveFactors(1,:),additiveFactors(3,:),additiveRange(k)); % slows down rate of radius growth appropriately
             
            Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
            Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

            % memory efficient high resolution solution
            [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
            averageRadius = m31./m21;
            averageHeight = m22./m21;
            stdRadius = sqrt(abs(m41./m21 - averageRadius.^2));
            stdHeight = sqrt(abs(m23./m21 - averageHeight.^2));
            temperatureRamp = [t; temperature]; % to be used by mom
            % method of moments solution:
            options = odeset('RelTol',1e-8,'AbsTol',1e-10);
            [t_mom, y] = ode113(@(t_mom, y)mom2D_additiveTD(t_mom,y,kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,temperatureRamp,shapeFactor,particleDensity,solubilityFactor,growthFactor), [0 simulationTime], y0,options);
            concentration_mom = y(:,17);
            m00_mom = y(:,1);
            m21_mom = y(:,7);
            m31_mom = y(:,11);
            m22_mom = y(:,9);
            m41_mom = y(:,16);
            m23_mom = y(:,15);
            averageRadius_mom = m31_mom./m21_mom;
            averageHeight_mom = m22_mom./m21_mom;
            stdRadius_mom = sqrt(abs(m41_mom./m21_mom - averageRadius_mom.^2));
            stdHeight_mom = sqrt(abs(m23_mom./m21_mom - averageHeight_mom.^2));

% aspect ratio plots:
            figure(1)
            if k == 1
                tiles1 = tiledlayout(gcf,1,3,'Units','inches','TileSpacing','compact','Padding','compact');
            end
            nexttile
            plot(averageRadius, averageHeight,"LineWidth",1)
            hold on
            
            % iso-volume plots:
            l1 = 100:0.1:1000;
            avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
            l2 = avgVolumemin./l1.^2;
            plot(l1,l2,'--',"LineWidth",1) % lower bound
            l2 = yieldFactor*avgVolumemin./l1.^2;
            plot(l1,l2,'--',"LineWidth",1) % desired final volume

            % ideal aspect ratio:
            plot(L2,L2,':',"LineWidth",1)
            xlim([200 600]);
            ylim([200 800]);

            % PSSD plots:
            contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.6); % initial PSSD
            contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.6); % final PSSD

            set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
            title(subplotTitles(k),'FontWeight','bold','FontSize',12)
            colororder(['#000000';'#d00000';'#03045e'; '#000000'])
            xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',10)
            xticks(linspace(200,600,3))
            ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',10)
            yticks(linspace(200,800,4))
            % pbaspect([1 1.1 1])
            

% concentration plot
            figure(2)
            plot(t,concentration,"LineWidth",1)
            hold on

% solubility plot
            figure(3)
            if k == 1
                tiles2 = tiledlayout(gcf,1,3,'Padding','compact','TileSpacing','compact','Units','inches');  
            end
            nexttile
            plot(temperature,concentration,"LineWidth",1)
            hold on
            plot(temperature(supersaturation>=1),concentration(supersaturation>=1),"LineWidth",1)
            
            % generate solubility curve for reference
            T = 0:0.1:40;
            solubility =solubilityFactor*3.37*exp(0.036*T);
            plot(T,solubility,"LineWidth",1)
            
            set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
            title(subplotTitles(k),'FontWeight','bold','Fontsize',12)
            xlabel({['{\it T} [' char(176) 'C]']},'FontSize',10)
            ylabel({'{\it c} [g kg^{-1}]'},'FontSize',10)
            colororder(gca,['#D00000';'#023E8A';'#000000'])
            xline([Tmin Tmax],'--',{'{\it T}_{min}','{\it T}_{max}'})
            % pbaspect([1 1 1])
            

% PSSD plots
            figure((k-1)*6+4)
            contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9]) % initial PSSD
            hold on
            contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'--') % final PSSD
            xlim([150 850]);
            ylim([150 850]);
    
            set(gca,'FontSize',8,'FontName','times')
            fig = gcf;
            fig.Units = "inches";
            fig.OuterPosition(3)=3.25;
            title('PSSD')
            xlabel({'L1', ['[' char(181) 'm]']})
            ylabel({'L2', ['[' char(181) 'm]']})

% standard deviation plots
            figure((k-1)*6+5)
            subplot(2,1,1)
            plot(t,stdRadius,"LineWidth",1)
            hold on
            plot(t_mom,stdRadius_mom,'--',"LineWidth",1)
    
            set(gca,'FontSize',8,'FontName','times')
            fig = gcf;
            fig.Units = "inches";
            fig.OuterPosition(3)=3.25;
            xlabel({'Time', '[h]'})
            ylabel({'\sigma_{11}', ['[' char(181) 'm]']})
            legend('High resolution', 'Method of moments')
    
            subplot(2,1,2)
            plot(t,stdHeight,"LineWidth",1)
            hold on
            plot(t_mom,stdHeight_mom,'--',"LineWidth",1)
    
            set(gca,'FontSize',8,'FontName','times')
            fig = gcf;
            fig.Units = "inches";
            fig.OuterPosition(3)=3.25;
            xlabel({'Time', '[h]'})
            ylabel({'\sigma_{22}', ['[' char(181) 'm]']})
            legend('High resolution', 'Method of moments')

% moment plots:
            figure(6)
            subplot(3,2,1)
            plot(t,m00,"LineWidth",1)
            hold on
            plot(t_mom,m00_mom,'--',"LineWidth",1)
    
            set(gca,'FontSize',20)
            title('m00 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,2)
            plot(t,m21,"LineWidth",1)
            hold on
            plot(t_mom,m21_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m21 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,3)
            plot(t,m31,"LineWidth",1)
            hold on
            plot(t_mom,m31_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m31 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,4)
            plot(t,m22,"LineWidth",1)
            hold on
            plot(t_mom,m22_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m22 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,5)
            plot(t,m41,"LineWidth",1)
            hold on
            plot(t_mom,m41_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m41 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,6)
            plot(t,m23,"LineWidth",1)
            hold on
            plot(t_mom,m23_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m23 profile')
            legend('High resolution','Method of moments')
        end

% add legends etc. to combined plots then save.
        figure(1)
        leg = legend('Average volume','Minimum average volume','Desired average volume','Ideal aspect ratio','Orientation','horizontal','FontSize',10,'Box','off');
        leg.Layout.Tile = 'south';

        tiles1.InnerPosition(1:2) = [3 3];
        tiles1.InnerPosition(4) = 2.25;
        tiles1.OuterPosition(3) = 7;

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];
          
        fileName = strcat('Aspect ratio plot for a) cA = ',num2str(additiveRange(1)),', b) cA = ',num2str(additiveRange(2)),', c) cA = ',num2str(additiveRange(3)),' [kg kg^-1]','.pdf');
        fpath = strcat(path,'Batch2D_additive\Additive variation cycles\', folderName);
        exportgraphics(tiles1, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

        figure(2)
        set(gca,'FontSize',8,'FontName','times','Units','inches')
        ax = gca;
        ax.InnerPosition(1:2) = [3 3];
        ax.InnerPosition(4) = 2;
        ax.OuterPosition(3) = 7;
        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];
        xlabel({'{\it t} [h]'},'Fontsize',10)
        ylabel({'{\it c} [g kg^{-1}]'},'Fontsize',10)
        legend(['{\it c}_A = ' num2str(additiveRange(1)) ' [kg kg^{-1}]'], ['{\it c}_A = ' num2str(additiveRange(2)) ' [kg kg^{-1}]'], ['{\it c}_A = ' num2str(additiveRange(3)) ' [kg kg^{-1}]'],'Box','off','Orientation','horizontal','FontSize',10,'Location','southoutside');
        colororder(['#FFBA08';'#DC2F02';'#6A040F']);
        h = get(gca,'Children');
        h(2).LineStyle = '--';
        h(3).LineStyle = ':';
        fileName = strcat('Concentration profiles for cA = ',num2str(additiveRange(1)),{' '},num2str(additiveRange(2)),{' '},num2str(additiveRange(3)),' [kg kg^-1]','.pdf');
        fpath = strcat(path,'Batch2D_additive\Additive variation cycles\', folderName);
        exportgraphics(gca, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

        figure(3)
        leg = legend('Dissolution stage','Growth stage', ...
                'Solubility curve','AutoUpdate','off','Orientation','horizontal','FontSize',10,'Box','off');
        leg.Layout.Tile = 'south';

        tiles2.InnerPosition(1:2) = [3 3];
        tiles2.InnerPosition(4) = 2.25;
        tiles2.OuterPosition(3) = 7;

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('Solubility plots for a) cA = ',num2str(additiveRange(1)),', b) cA = ',num2str(additiveRange(2)),', c) cA = ',num2str(additiveRange(3)),' [kg kg^-1]','.pdf');
        fpath = strcat(path,'Batch2D_additive\Additive variation cycles\', folderName);
        exportgraphics(tiles2, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

% store operating conditions etc in a matrix to be exported as a
% text file
        namesCell = {'Initial Concentration';'Shape factor'; 'Particle density'; 'Simulation time';'dL1';'dL2';'Yield';'Number of cycles';'Minimum supersaturation';'Maximum supersaturation'};
        unitsCell = {'[g/kg]'; '[-]';['[g/' char(181) 'm' char(179) ']'];'[h]';['[' char(181) 'm]'];['[' char(181) 'm]'];'[-]';'[-]';'[-]';'[-]'};
        parametersCell = cat(1, [namesCell, [essentialParameters(1); essentialParameters(3:7); cycleParameters], unitsCell; {'%', 'Additive concentration range:','%'}; num2cell(additiveRange)]);
        fileName='Additive variation parameters';
        fpath = strcat(path, 'Batch2D_additive\Additive variation cycles\', folderName);
        writecell(parametersCell,fullfile(fpath,fileName))
%% Supersaturation effect on temperature cycles
    case 3
        % create a folder with the current date and time:
        mkdir(strcat(path,'Batch2D_additive\Supersaturation variation cycles'), folderName)

% prompt user for desired volume and supersaturation limits
        prompt = {'Desired final yield (as a multiple of the original mass):','Number of growth/dissolution cycles:','Minimum supersaturation:','Maximum supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'3','2','0.9','1.15'};
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
        Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
        Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

        supersaturationRange = [1.1 1.15 1.2];
        subplotTitles = {'a)', 'b)', 'c)'};
        for k = 1:length(supersaturationRange)
            Smax = supersaturationRange(k);
            supersaturationLimits = [Smin Smax];

% memory efficient high resolution solution
            [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
            averageRadius = m31./m21;
            averageHeight = m22./m21;
            stdRadius = sqrt(abs(m41./m21 - averageRadius.^2));
            stdHeight = sqrt(abs(m23./m21 - averageHeight.^2));
            temperatureRamp = [t; temperature]; % to be used by mom
% method of moments solution:
            options = odeset('RelTol',1e-8,'AbsTol',1e-10);
            [t_mom, y] = ode113(@(t_mom, y)mom2D_additiveTD(t_mom,y,kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,temperatureRamp,shapeFactor,particleDensity,solubilityFactor,growthFactor), [0 simulationTime], y0,options);
            concentration_mom = y(:,17);
            m00_mom = y(:,1);
            m21_mom = y(:,7);
            m31_mom = y(:,11);
            m22_mom = y(:,9);
            m41_mom = y(:,16);
            m23_mom = y(:,15);
            averageRadius_mom = m31_mom./m21_mom;
            averageHeight_mom = m22_mom./m21_mom;
            stdRadius_mom = sqrt(abs(m41_mom./m21_mom - averageRadius_mom.^2));
            stdHeight_mom = sqrt(abs(m23_mom./m21_mom - averageHeight_mom.^2));

% aspect ratio plots:
            figure(1)
            if k == 1
                tiles1 = tiledlayout(1,3,'Units','inches','TileSpacing','compact','Padding','compact');
            end
            nexttile
            plot(averageRadius, averageHeight,"LineWidth",1)
            hold on
            
            % iso-volume plots:
            l1 = 100:0.1:1000;
            avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
            l2 = avgVolumemin./l1.^2;
            plot(l1,l2,'--',"LineWidth",1) % lower bound
            l2 = yieldFactor*avgVolumemin./l1.^2;
            plot(l1,l2,'--',"LineWidth",1) % desired final volume

            % ideal aspect ratio:
            plot(L2,L2,':',"LineWidth",1)
            xlim([200 600]);
            ylim([200 800]);

            % PSSD plots:
            contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.6); % initial PSSD
            contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.6); % final PSSD

            set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
            title(subplotTitles(k),'FontWeight','bold','FontSize',12)
            colororder(['#000000';'#d00000';'#03045e'; '#000000'])
            xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',10)
            xticks(linspace(200,600,3))
            ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',10)
            yticks(linspace(200,800,4))
            % pbaspect([1 1.1 1])
            

% concentration plot
            figure(2)
            plot(t,concentration,"LineWidth",1)
            hold on

% solubility plot
            figure(3)
            if k == 1
                tiles2 = tiledlayout(1,3,'Padding','compact','TileSpacing','compact','Units','inches');  
            end
            nexttile
            plot(temperature,concentration,"LineWidth",1)
            hold on
            plot(temperature(supersaturation>=1),concentration(supersaturation>=1),"LineWidth",1)
            
            % generate solubility curve for reference
            T = 0:0.1:40;
            solubility =solubilityFactor*3.37*exp(0.036*T);
            plot(T,solubility,"LineWidth",1)
            
            set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
            title(subplotTitles(k),'FontWeight','bold','Fontsize',12)
            xlabel({['{\it T} [' char(176) 'C]']},'FontSize',10)
            ylabel({'{\it c} [g kg^{-1}]'},'FontSize',10)
            colororder(['#D00000';'#023E8A';'#000000'])
            xline([Tmin Tmax],'--',{'{\it T}_{min}','{\it T}_{max}'})
            % pbaspect([1 1 1])


% PSSD plots
            figure((k-1)*6+4)
            contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9]) % initial PSSD
            hold on
            contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'--') % final PSSD
            xlim([150 850]);
            ylim([150 850]);
    
            set(gca,'FontSize',8,'FontName','times')
            fig = gcf;
            fig.Units = "inches";
            fig.OuterPosition(3)=3.25;
            title('PSSD')
            xlabel({'L1', ['[' char(181) 'm]']})
            ylabel({'L2', ['[' char(181) 'm]']})

% standard deviation plots
            figure((k-1)*6+5)
            subplot(2,1,1)
            plot(t,stdRadius,"LineWidth",1)
            hold on
            plot(t_mom,stdRadius_mom,'--',"LineWidth",1)
    
            set(gca,'FontSize',8,'FontName','times')
            fig = gcf;
            fig.Units = "inches";
            fig.OuterPosition(3)=3.25;
            xlabel({'Time', '[h]'})
            ylabel({'\sigma_{11}', ['[' char(181) 'm]']})
            legend('High resolution', 'Method of moments')
    
            subplot(2,1,2)
            plot(t,stdHeight,"LineWidth",1)
            hold on
            plot(t_mom,stdHeight_mom,'--',"LineWidth",1)
    
            set(gca,'FontSize',8,'FontName','times')
            fig = gcf;
            fig.Units = "inches";
            fig.OuterPosition(3)=3.25;
            xlabel({'Time', '[h]'})
            ylabel({'\sigma_{22}', ['[' char(181) 'm]']})
            legend('High resolution', 'Method of moments')

% moment plots:
            figure(6)
            subplot(3,2,1)
            plot(t,m00,"LineWidth",1)
            hold on
            plot(t_mom,m00_mom,'--',"LineWidth",1)
    
            set(gca,'FontSize',20)
            title('m00 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,2)
            plot(t,m21,"LineWidth",1)
            hold on
            plot(t_mom,m21_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m21 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,3)
            plot(t,m31,"LineWidth",1)
            hold on
            plot(t_mom,m31_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m31 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,4)
            plot(t,m22,"LineWidth",1)
            hold on
            plot(t_mom,m22_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m22 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,5)
            plot(t,m41,"LineWidth",1)
            hold on
            plot(t_mom,m41_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m41 profile')
            legend('High resolution','Method of moments')
    
            subplot(3,2,6)
            plot(t,m23,"LineWidth",1)
            hold on
            plot(t_mom,m23_mom,'--',"LineWidth",1)
            set(gca,'FontSize',20)
            title('m23 profile')
            legend('High resolution','Method of moments')
        end

% add legends etc. to combined plots then save.
        figure(1)
        leg = legend('Average volume','Minimum average volume','Desired average volume','Ideal aspect ratio','Orientation','horizontal','FontSize',10,'Box','off');
        leg.Layout.Tile = 'south';

        tiles1.InnerPosition(1:2) = [3 3];
        tiles1.InnerPosition(4) = 2.25;
        tiles1.OuterPosition(3) = 7;

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('Aspect ratio plots for a) S = ',num2str(supersaturationRange(1)),', b) S = ',num2str(supersaturationRange(2)), ' c) S = ',num2str(supersaturationRange(3)), '[-]','.pdf');
        fpath = strcat(path,'Batch2D_additive\Supersaturation variation cycles\', folderName);
        exportgraphics(tiles1, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

        figure(2)
        set(gca,'FontSize',8,'FontName','times','Units','inches')
        ax = gca;
        ax.InnerPosition(1:2) = [3 3];
        ax.InnerPosition(4) = 2;
        ax.OuterPosition(3) = 7;
        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];
        xlabel({'{\it t} [h]'},'Fontsize',10)
        ylabel({'{\it c} [g kg^{-1}]'},'Fontsize',10)
        legend(['{\it S} = ' num2str(supersaturationRange(1))], ['{\it S} = ' num2str(supersaturationRange(2))], ['{\it S} = ' num2str(supersaturationRange(3))],'Box','off','Orientation','horizontal','FontSize',10,'Location','southoutside')
        colororder(['#FFBA08';'#DC2F02';'#6A040F']);
        h = get(gca,'Children');
        h(2).LineStyle = '--';
        h(3).LineStyle = ':';
        fileName = strcat('Concentration profiles for S = ',num2str(supersaturationRange(1)),{' '},num2str(supersaturationRange(2)),{' '},num2str(supersaturationRange(3)),' [-]','.pdf');
        fpath = strcat(path,'Batch2D_additive\Supersaturation variation cycles\', folderName);
        exportgraphics(gca, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

        figure(3)
        leg = legend('Dissolution stage','Growth stage', ...
                'Solubility curve','AutoUpdate','off','Orientation','horizontal','FontSize',10,'Box','off');
        leg.Layout.Tile = 'south';

        tiles2.InnerPosition(1:2) = [3 3];
        tiles2.InnerPosition(4) = 2.25;
        tiles2.OuterPosition(3) = 7;

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('Solubility plots for a) S = ',num2str(supersaturationRange(1)),', b) S = ',num2str(supersaturationRange(2)), ' c) S = ',num2str(supersaturationRange(3)), '[-]','.pdf');
        fpath = strcat(path,'Batch2D_additive\Supersaturation variation cycles\', folderName);
        exportgraphics(tiles2, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');
        
% store operating conditions etc in a matrix to be exported as a
% text file
        namesCell = {'Initial Concentration';'additive concentration [kg/kg]s';'Shape factor'; 'Particle density'; 'Simulation time';'dL1';'dL2';'Yield';'Number of cycles';'Minimum supersaturation'};
        unitsCell = {'[g/kg]'; '[-]';['[g/' char(181) 'm' char(179) ']'];'[h]';['[' char(181) 'm]'];['[' char(181) 'm]'];'[-]';'[-]';'[-]';'[-]'};
        parametersCell = cat(1, [namesCell, [essentialParameters; cycleParameters(1:3)], unitsCell; {'%', 'Supersaturation range:','%'}; num2cell(supersaturationRange)]);
        fileName='Supersaturation variation parameters';
        fpath = strcat(path,'Batch2D_additive\Supersaturation variation cycles\', folderName);
        writecell(parametersCell,fullfile(fpath,fileName))

%% Relative kinetics variation
    case 4
% create a folder with the current date and time:
        mkdir(strcat(path,'Batch2D_additive\Relative kinetics variation cycles'), folderName)

% prompt user for desired volume and supersaturation limits
        prompt = {'Desired final yield (as a multiple of the original mass):','Number of growth/dissolution cycles:','Minimum supersaturation:','Maximum supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'3','2','0.9','1.15'};
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
        Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
        Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

        growthRange = [0.5 1 2];
        dissolutionRange = [0.75 1 4];
        subplotTitles = {'a)', 'b)', 'c)'};
        for i = 1:length(dissolutionRange)
            for k = 1:length(growthRange)
                kg21var = growthRange(k)*kg21;
                kd21var = dissolutionRange(i)*kd21;
    
% memory efficient high resolution solution
                [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21var, kg22, kg23, kd11, kd12, kd21var, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
                averageRadius = m31./m21;
                averageHeight = m22./m21;
                stdRadius = sqrt(abs(m41./m21 - averageRadius.^2));
                stdHeight = sqrt(abs(m23./m21 - averageHeight.^2));
                temperatureRamp = [t; temperature]; % to be used by mom
% method of moments solution:
                options = odeset('RelTol',1e-8,'AbsTol',1e-10);
                [t_mom, y] = ode113(@(t_mom, y)mom2D_additiveTD(t_mom,y,kg11,kg12,kg13,kg21var,kg22,kg23,kd11,kd12,kd21var,kd22,temperatureRamp,shapeFactor,particleDensity,solubilityFactor,growthFactor), [0 simulationTime], y0,options);
                concentration_mom = y(:,17);
                m00_mom = y(:,1);
                m21_mom = y(:,7);
                m31_mom = y(:,11);
                m22_mom = y(:,9);
                m41_mom = y(:,16);
                m23_mom = y(:,15);
                averageRadius_mom = m31_mom./m21_mom;
                averageHeight_mom = m22_mom./m21_mom;
                stdRadius_mom = sqrt(abs(m41_mom./m21_mom - averageRadius_mom.^2));
                stdHeight_mom = sqrt(abs(m23_mom./m21_mom - averageHeight_mom.^2));
    
% aspect ratio plots:
                figure(1)
                if k == 1 && i == 1
                    tiles1 = tiledlayout(gcf,3,3,'Units','inches','TileSpacing','compact','Padding','compact');
                end
                nexttile
                plot(averageRadius, averageHeight,"LineWidth",1)
                hold on
                
                % iso-volume plots:
                l1 = 100:0.1:1000;
                avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
                l2 = avgVolumemin./l1.^2;
                plot(l1,l2,'--',"LineWidth",1) % lower bound
                l2 = yieldFactor*avgVolumemin./l1.^2;
                plot(l1,l2,'--',"LineWidth",1) % desired final volume
    
                % ideal aspect ratio:
                plot(L2,L2,':',"LineWidth",1)
                xlim([200 600]);
                ylim([200 900]);
    
                % PSSD plots:
                contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.6); % initial PSSD
                contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.6); % final PSSD
    
                set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
                colororder(['#000000';'#d00000';'#03045e'; '#000000'])
                xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',10)
                xticks(linspace(200,600,3))
                ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',10)
                yticks(linspace(200,900,5))
                % pbaspect([1 1 1])
                
    
% concentration plot
                figure(2)

                if k == 1 && i == 1
                    tiles3 = tiledlayout(gcf,2,1,'Padding','compact','TileSpacing','compact','Units','inches');  
                end

                if i == 2 % G2 concentration plot  

                    nexttile(1)

                    plot(t,concentration,"LineWidth",1)
                    hold on

                    set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
                    colororder(gca,['#FFBA08';'#E85D04';'#D00000']);
                    h = get(gca,'Children');
                    if length(h) == 3
                        h(2).LineStyle = '--';
                        h(3).LineStyle = ':';
                    end
                    xlabel({'{\it t} [h]'},'Fontsize',10)
                    ylabel({'{\it c} [g kg^{-1}]'},'Fontsize',10)
                end

                if k == 2 % D2 concentration plot
                    nexttile(2)

                    plot(t,concentration,"LineWidth",1)
                    hold on

                    set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
                    colororder(gca,['#00B4D8';'#0077B6';'#03045E']);
                    h = get(gca,'Children');
                    if length(h) == 3
                        h(2).LineStyle = '--';
                        h(3).LineStyle = ':';
                    end
                    xlabel({'{\it t} [h]'},'Fontsize',10)
                    ylabel({'{\it c} [g kg^{-1}]'},'Fontsize',10)

                end

% solubility plot
                figure(3)
                if k == 1 && i == 1
                    tiles2 = tiledlayout(gcf,3,3,'Padding','compact','TileSpacing','compact','Units','inches');  
                end

                nexttile
                plot(temperature,concentration,"LineWidth",1)
                hold on
                plot(temperature(supersaturation<1),concentration(supersaturation<1),"LineWidth",1)
                
                % generate solubility curve for reference
                T = 0:0.1:40;
                solubility =solubilityFactor*3.37*exp(0.036*T);
                plot(T,solubility,"LineWidth",1)
                
                h = get(gca,'Children');
                set(gca,'FontSize',8,'FontName','times','Children',[h(1) h(2) h(3)],'TitleHorizontalAlignment','left')
                xlabel({['{\it T} [' char(176) 'C]']},'FontSize',10)
                ylabel({'{\it c} [g kg^{-1}]'},'FontSize',10)
                colororder(['#023E8A';'#D00000';'#000000'])
                xline([Tmin Tmax],'--',{'{\it T}_{min}','{\it T}_{max}'})
                % pbaspect([1 1 1])
    
    
% PSSD plots
                figure((k-1)*6+4)
                contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9]) % initial PSSD
                hold on
                contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'--') % final PSSD
                xlim([150 850]);
                ylim([150 850]);
        
                set(gca,'FontSize',8,'FontName','times')
                fig = gcf;
                fig.Units = "inches";
                fig.OuterPosition(3)=3.25;
                title('PSSD')
                xlabel({'L1', ['[' char(181) 'm]']})
                ylabel({'L2', ['[' char(181) 'm]']})
    
% standard deviation plots
                figure((k-1)*6+5)
                subplot(2,1,1)
                plot(t,stdRadius,"LineWidth",1)
                hold on
                plot(t_mom,stdRadius_mom,'--',"LineWidth",1)
        
                set(gca,'FontSize',8,'FontName','times')
                fig = gcf;
                fig.Units = "inches";
                fig.OuterPosition(3)=3.25;
                xlabel({'Time', '[h]'})
                ylabel({'\sigma_{11}', ['[' char(181) 'm]']})
                legend('High resolution', 'Method of moments')
        
                subplot(2,1,2)
                plot(t,stdHeight,"LineWidth",1)
                hold on
                plot(t_mom,stdHeight_mom,'--',"LineWidth",1)
        
                set(gca,'FontSize',8,'FontName','times')
                fig = gcf;
                fig.Units = "inches";
                fig.OuterPosition(3)=3.25;
                xlabel({'Time', '[h]'})
                ylabel({'\sigma_{22}', ['[' char(181) 'm]']})
                legend('High resolution', 'Method of moments')
    
% moment plots:
                figure(6)
                subplot(3,2,1)
                plot(t,m00,"LineWidth",1)
                hold on
                plot(t_mom,m00_mom,'--',"LineWidth",1)
        
                set(gca,'FontSize',20)
                title('m00 profile')
                legend('High resolution','Method of moments')
        
                subplot(3,2,2)
                plot(t,m21,"LineWidth",1)
                hold on
                plot(t_mom,m21_mom,'--',"LineWidth",1)
                set(gca,'FontSize',20)
                title('m21 profile')
                legend('High resolution','Method of moments')
        
                subplot(3,2,3)
                plot(t,m31,"LineWidth",1)
                hold on
                plot(t_mom,m31_mom,'--',"LineWidth",1)
                set(gca,'FontSize',20)
                title('m31 profile')
                legend('High resolution','Method of moments')
        
                subplot(3,2,4)
                plot(t,m22,"LineWidth",1)
                hold on
                plot(t_mom,m22_mom,'--',"LineWidth",1)
                set(gca,'FontSize',20)
                title('m22 profile')
                legend('High resolution','Method of moments')
        
                subplot(3,2,5)
                plot(t,m41,"LineWidth",1)
                hold on
                plot(t_mom,m41_mom,'--',"LineWidth",1)
                set(gca,'FontSize',20)
                title('m41 profile')
                legend('High resolution','Method of moments')
        
                subplot(3,2,6)
                plot(t,m23,"LineWidth",1)
                hold on
                plot(t_mom,m23_mom,'--',"LineWidth",1)
                set(gca,'FontSize',20)
                title('m23 profile')
                legend('High resolution','Method of moments')
            end
        end

% add legends etc. to combined plots then save.
        figure(1)
        leg = legend('Average volume','Minimum average volume','Desired average volume','Ideal aspect ratio','Orientation','horizontal','FontSize',10,'Box','off');
        leg.Layout.Tile = 'south';

        tiles1.InnerPosition(1:2) = [3 3];
        tiles1.InnerPosition(4) = 6.25;
        tiles1.OuterPosition(3) = 7;

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('Aspect ratio plots for increasing G2 (left to right) and increasing D2 (top to bottom).pdf');
        fpath = strcat(path,'Batch2D_additive\Relative kinetics variation cycles\', folderName);
        exportgraphics(tiles1, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

        figure(2)

        tiles3.InnerPosition(1:2) = [3 3];
        tiles3.InnerPosition(4) = 4;
        tiles3.OuterPosition(3) = 7;
        
        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];
        
            
        fileName = strcat('Concentration profiles for increasing G2 (top) and increasing D2 (bottom).pdf');
        fpath = strcat(path,'Batch2D_additive\Relative kinetics variation cycles\', folderName);
        exportgraphics(tiles3, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

        figure(3)
        leg = legend('Growth stage','Dissolution stage', ...
                'Solubility curve','AutoUpdate','off','Orientation','horizontal','FontSize',10,'Box','off');
        leg.Layout.Tile = 'south';

        tiles2.InnerPosition(1:2) = [3 3];
        tiles2.InnerPosition(4) = 6.25;
        tiles2.OuterPosition(3) = 7;

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('Solubility plots for increasing G2 (left to right) and increasing D2 (top to bottom).pdf');
        fpath = strcat(path,'Batch2D_additive\Relative kinetics variation cycles\', folderName);
        exportgraphics(tiles2, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');
        
% store operating conditions etc in a matrix to be exported as a
% text file
        namesCell = {'Initial Concentration';'additive concentration [kg/kg]s';'Shape factor'; 'Particle density'; 'Simulation time';'dL1';'dL2';'Yield';'Number of cycles';'Minimum supersaturation';'Maximum supersaturation'};
        unitsCell = {'[g/kg]'; '[-]';['[g/' char(181) 'm' char(179) ']'];'[h]';['[' char(181) 'm]'];['[' char(181) 'm]'];'[-]';'[-]';'[-]';'[-]';'[-]'};
        parametersCell = cat(1, [namesCell, [essentialParameters; cycleParameters], unitsCell; {'%', 'Growth range:','%'}; num2cell(growthRange); {'%', 'Dissolution range:','%'}; num2cell(dissolutionRange)]);
        fileName='Relative kinetics variation parameters';
        fpath = strcat(path,'Batch2D_additive\Relative kinetics variation cycles\', folderName);
        writecell(parametersCell,fullfile(fpath,fileName))

%% Pure growth base case
    case 5
% create a folder with the current date and time:
        mkdir(strcat(path,'Batch2D_additive\Pure growth base case'), folderName)

% prompt user for desired volume and supersaturation limits
        prompt = {'Desired final yield (as a multiple of the original mass):','Constant supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'3','1.15'};
        cycleParameters = inputdlg(prompt,dlgtitle, [1 35],definput);
        yieldFactor = str2double(cycleParameters(1));
        maxCycleNumber = 0;
        Smin = 0.9;
        Smax = str2double(cycleParameters(2));
        supersaturationLimits = [Smin Smax];
        m21min = initialm21;
        m21max = yieldFactor*m21min;
        cmax = initialConcentration;
        cmin = cmax - particleDensity*shapeFactor*(m21max - m21min);
        Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
        Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

% memory efficient high resolution solution:
        [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
        averageRadius = m31./m21;
        averageHeight = m22./m21;
        stdRadius = sqrt(abs(m41./m21 - averageRadius.^2));
        stdHeight = sqrt(abs(m23./m21 - averageHeight.^2));
        temperatureRamp = [t; temperature]; % to be used by mom
        
% method of moments solution:
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t_mom, y] = ode113(@(t_mom, y)mom2D_additiveTD(t_mom,y,kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,temperatureRamp,shapeFactor,particleDensity,solubilityFactor,growthFactor), [0 simulationTime], y0,options);
        concentration_mom = y(:,17);
        m00_mom = y(:,1);
        m21_mom = y(:,7);
        m31_mom = y(:,11);
        m22_mom = y(:,9);
        m41_mom = y(:,16);
        m23_mom = y(:,15);
        averageRadius_mom = m31_mom./m21_mom;
        averageHeight_mom = m22_mom./m21_mom;
        stdRadius_mom = sqrt(abs(m41_mom./m21_mom - averageRadius_mom.^2));
        stdHeight_mom = sqrt(abs(m23_mom./m21_mom - averageHeight_mom.^2));

% aspect ratio plots:
        figure(1)
        tiles1 = tiledlayout(1,2,'Padding','compact','TileSpacing','compact','Units','inches');
        tiles1.InnerPosition(1:2) = [3 3];
        tiles1.InnerPosition(4) = 3.25;
        tiles1.OuterPosition(3) = 7;
        nexttile
        plot(averageRadius, averageHeight,"LineWidth",1)
        hold on
        plot(averageRadius_mom, averageHeight_mom, '--',"LineWidth",1)
        
        % iso-volume plots:
        l1 = 100:0.1:1000;
        avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
        l2 = avgVolumemin./l1.^2;
        plot(l1,l2,"LineWidth",1) % lower bound
        l2 = yieldFactor*avgVolumemin./l1.^2;
        plot(l1,l2,"LineWidth",1) % desired final volume

        % ideal aspect ratio:
        plot(L2,L2,':',"LineWidth",1)
        xlim([200 600]);
        ylim([200 800]);

        % PSSD plots
        contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.6); % initial PSSD
        contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.6); % final PSSD

        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        colororder(gca,['#000000'; '#00A300';'#d00000';'#03045e'; '#000000'])
        xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',10)
        xticks(linspace(200,600,3))
        ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',10)
        yticks(linspace(200,800,4))
        % pbaspect([1 1.1 1])
        legend('High resolution solution','Method of moments solution','Minimum volume','Desired volume','Ideal aspect ratio','Location','southeast')
        title('a)','FontWeight','bold','Fontsize',12)
        text([averageRadius(1) averageRadius(end)], [averageHeight(1) averageHeight(end)],{'     start' '     end'})
        
% solubility plot:
        figure(1)
        nexttile
        plot(temperature,concentration,"LineWidth",1)
        hold on

        % generate solubility curve for reference
        T = 0:0.1:40;
        solubility =solubilityFactor*3.37*exp(0.036*T);
        plot(T,solubility,"LineWidth",1)
        
        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        xlabel({['{\it T} [' char(176) 'C]']},'FontSize',10)
        ylabel({'{\it c} [g kg^{-1}]'},'FontSize',10)
        % pbaspect([1 1.1 1])
        legend('Growth stage', ...
            'Solubility curve','AutoUpdate','off','Location','southeast');
        xline([Tmin Tmax],'--',{'{\it T}_{min}','{\it T}_{max}'})
        title('b)','FontWeight','bold','Fontsize',12)
        colororder(gca,['#023E8A';'#000000'])
        text([temperature(1) temperature(end)], [concentration(1) concentration(end)],{'\leftarrow start' '\leftarrow end'})

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('aspect ratio + solubility base case.pdf');
        fpath = strcat(path,'Batch2D_additive\Pure growth base case\', folderName);
        exportgraphics(tiles1, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

% concentration plot:
        figure(2)
        tiles2 = tiledlayout(3,1,'Padding','compact','TileSpacing','compact','Units','inches');
        tiles2.InnerPosition(1:2) = [3 3];
        tiles2.InnerPosition(4) = 6;
        tiles2.OuterPosition(3) = 7;
        nexttile
        plot(t,concentration,"LineWidth",1)
        hold on
        plot(t_mom,concentration_mom,'--',"LineWidth",1)

        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        xlabel({'{\it t} [h]'},'FontSize',10)
        ylabel({'{\it c} [g kg^{-1}]'},'FontSize',10)
        title('a)','FontWeight','bold','Fontsize',12)
        colororder(gca,['#000000';'#00A300'])
        leg = legend('High resolution solution', 'Method of moments solution','AutoUpdate','off','Orientation','horizontal','FontSize',10,'Box','off');
% standard deviation plots
        nexttile
        plot(t,stdRadius./stdRadius(1),"LineWidth",1)
        hold on
        plot(t_mom,stdRadius_mom./stdRadius_mom(1),'--',"LineWidth",1)

        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        xlabel({'{\it t} [h]'},'FontSize',10)
        ylabel({'\sigma_{11,model}/\sigma_{11,actual} [-]'},'FontSize',10)
        title('b)','FontWeight','bold','Fontsize',12)
        colororder(gca,['#000000';'#00A300'])

        nexttile
        plot(t,stdHeight./stdHeight(1),"LineWidth",1)
        hold on
        plot(t_mom,stdHeight_mom./stdHeight_mom(1),'--',"LineWidth",1)

        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
        xlabel({'{\it t} [h]'},'FontSize',10)
        ylabel({'\sigma_{22,model}/\sigma_{22,actual} [-]'},'FontSize',10)
        title('c)','FontWeight','bold','Fontsize',12)
        colororder(gca,['#000000';'#00A300'])
        
        leg.Layout.Tile = 'south';

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];
        
        fileName = strcat('Concentration + standard deviation profiles base case.pdf');
        fpath = strcat(path,'Batch2D_additive\Pure growth base case\', folderName);
        exportgraphics(tiles2, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

% PSSD plots
        figure(3)
        contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9]) % initial PSSD
        hold on
        contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'--') % final PSSD
        xlim([150 850]);
        ylim([150 850]);

        set(gca,'FontSize',8,'FontName','times')
        fig = gcf;
        fig.Units = "inches";
        fig.OuterPosition(3)=3.25;
        title('PSSD')
        xlabel({'L1', ['[' char(181) 'm]']})
        ylabel({'L2', ['[' char(181) 'm]']})

% moment plots:
        figure(4)
        subplot(3,2,1)
        plot(t,m00,"LineWidth",1)
        hold on
        plot(t_mom,m00_mom,'--',"LineWidth",1)

        set(gca,'FontSize',20)
        title('m00 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,2)
        plot(t,m21,"LineWidth",1)
        hold on
        plot(t_mom,m21_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m21 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,3)
        plot(t,m31,"LineWidth",1)
        hold on
        plot(t_mom,m31_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m31 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,4)
        plot(t,m22,"LineWidth",1)
        hold on
        plot(t_mom,m22_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m22 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,5)
        plot(t,m41,"LineWidth",1)
        hold on
        plot(t_mom,m41_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m41 profile')
        legend('High resolution','Method of moments')

        subplot(3,2,6)
        plot(t,m23,"LineWidth",1)
        hold on
        plot(t_mom,m23_mom,'--',"LineWidth",1)
        set(gca,'FontSize',20)
        title('m23 profile')
        legend('High resolution','Method of moments')

% store operating conditions etc in a matrix to be exported as a
% text file
        namesCell = {'Initial Concentration'; 'Additive concentration';'Shape factor'; 'Particle density'; 'Simulation time';'dL1';'dL2';'Yield';'Constant supersaturation'};
        unitsCell = {'[g/kg]';'[kg/kg]'; '[-]';['[g/' char(181) 'm' char(179) ']'];'[h]';['[' char(181) 'm]'];['[' char(181) 'm]'];'[-]';'[-]'};
        parametersCell = cat(1, [namesCell, [essentialParameters(:); cycleParameters], unitsCell]);
        fileName='Pure growth base case parameters';
        fpath = strcat(path,'Batch2D_additive\Pure growth base case\', folderName);
        writecell(parametersCell,fullfile(fpath,fileName))

        %% Pure growth bin size variation
        case 6
% create a folder with the current date and time:
        mkdir(strcat(path,'Batch2D_additive\Bin size variation growth'), folderName)

% prompt user for desired volume and supersaturation limits
        prompt = {'Desired final yield (as a multiple of the original mass):','Constant supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'3','1.15'};
        cycleParameters = inputdlg(prompt,dlgtitle, [1 35],definput);
        yieldFactor = str2double(cycleParameters(1));
        maxCycleNumber = 0;
        Smin = 0.9;
        Smax = str2double(cycleParameters(2));
        supersaturationLimits = [Smin Smax];
        m21min = initialm21;
        m21max = yieldFactor*m21min;
        cmax = initialConcentration;
        cmin = cmax - particleDensity*shapeFactor*(m21max - m21min);
        Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
        Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

        % initialise
        startTime = zeros(3,1);
        endTime = zeros(3,1);
        stdRadius = cell(1,3);
        stdHeight = cell(1,3);
        t = cell(1,3);

        binRange = [1 2 5];
        for k = 1:length(binRange)
            dL1 = binRange(k);
            dL2 = binRange(k);

            L1 = 1:dL1:600; % [um]
            L2 = 1:dL2:1000; % [um]
            [l1, l2] = meshgrid(L1,L2);
            
            L = [l1(:) l2(:)];
            clear l1 l2
            meanDimensions = [250 400];
            standardDeviation = [50 0; 0 50];
            
            initialPSD = 1e5*mvnpdf(L,meanDimensions, standardDeviation);
            initialPSD = reshape(initialPSD,length(L2),length(L1));
             
            % memory efficient high resolution solution
            tic
            [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t{k}, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
            averageRadius = m31./m21;
            averageHeight = m22./m21;
            stdRadius{k} = sqrt(abs(m41./m21 - averageRadius.^2));
            stdHeight{k} = sqrt(abs(m23./m21 - averageHeight.^2));
            endTime(k) = toc;
        end

% standard deviation plots
        figure(1)

        tiles = tiledlayout(2,1,'Padding','compact','TileSpacing','compact','Units','inches');
        tiles.InnerPosition(1:2) = [3 3];
        tiles.InnerPosition(4) = 4;
        tiles.OuterPosition(3) = 7;
        nexttile
        plot(t{1},stdRadius{1}./stdRadius{1}(1),"LineWidth",1,'Color','#000000')
        hold on
        plot(t{2},stdRadius{2}./stdRadius{2}(1),"LineWidth",1,'Color','#000000')
        plot(t{3},stdRadius{3}./stdRadius{3}(1),"LineWidth",1,'Color','#000000')
        
        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left','LineStyleCyclingMethod','beforecolor')
        ax = gca;
        ax.LineStyleOrder = ["-";"--";":"];
        title('a)','FontWeight','bold','Fontsize',12)
        xlabel({'{\it t} [h]'},'FontSize',10)
        ylabel({'\sigma_{11,model}/\sigma_{11,actual} [-]'},'FontSize',10)

        nexttile
        plot(t{1},stdHeight{1}./stdHeight{1}(1),"LineWidth",1,'Color','#000000')
        hold on
        plot(t{2},stdHeight{2}./stdHeight{2}(1),"LineWidth",1,'Color','#000000')
        plot(t{3},stdHeight{3}./stdHeight{3}(1),"LineWidth",1,'Color','#000000')

        set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left','LineStyleCyclingMethod','beforecolor')
        ax = gca;
        ax.LineStyleOrder = ["-";"--";":"];
        title('b)','FontWeight','bold','Fontsize',12)
        xlabel({'{\it t} [h]'},'FontSize',10)
        ylabel({'\sigma_{22,model}/\sigma_{22,actual} [-]'},'FontSize',10)

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        leg = legend(['Bin size = ' num2str(binRange(1)) [' [' char(181) 'm]']], ['Bin size = ' num2str(binRange(2)) [' [' char(181) 'm]']], ['Bin size = ' num2str(binRange(3)) [' [' char(181) 'm]']],'Box','off','Orientation','horizontal','FontSize',10);
        leg.Layout.Tile = "south";

        fileName = strcat('Stdev profiles for bin size = ',num2str(binRange(1)),{' '},num2str(binRange(2)),{' '},num2str(binRange(3)),[' [' char(181) 'm]'],'.pdf');
        fpath = strcat(path,'Batch2D_additive\Bin size variation growth\', folderName);
        exportgraphics(tiles, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');


% store operating conditions etc in a matrix to be exported as a
% text file
        namesCell = {'Initial Concentration'; 'Additive concentration';'Shape factor'; 'Particle density'; 'Simulation time';'Yield';'Constant supersaturation'};
        unitsCell = {'[g/kg]';'[kg/kg]'; '[-]';['[g/' char(181) 'm' char(179) ']'];'[h]';'[-]';'[-]'};
        parametersCell = cat(1, [namesCell, [essentialParameters(1:5); cycleParameters], unitsCell; {'%', 'Bin size range [um]:','%'}; num2cell(binRange)]);
        fileName='Bin size variation parameters';
        fpath = strcat(path, 'Batch2D_additive\Bin size variation growth\', folderName);
        writecell(parametersCell,fullfile(fpath,fileName))

% calculate run times then print to cell
        namesCell = {'Time for bin size of 1 um in seconds = ';'Time for bin size of 2 um in seconds = ';'Time for bin size of 5 um in seconds = '};
        benchmark = cat(1, [namesCell, num2cell(endTime)]);
        fileName='Bin size variation benchmark';
        fpath = strcat(path, 'Batch2D_additive\Bin size variation growth\', folderName);
        writecell(benchmark,fullfile(fpath,fileName))

%% Attainable region pure growth
    case 8
        simulationTime = 2000; % [h]
% create a folder with the current date and time:
        mkdir(strcat(path,'Batch2D_additive\Attainable region growth'), folderName)
        
% prompt user for desired volume and supersaturation limits
        prompt = {'Desired final yield (as a multiple of the original mass):','Constant supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'3','1.15'};
        cycleParameters = inputdlg(prompt,dlgtitle, [1 35],definput);
        yieldFactor = str2double(cycleParameters(1));
        maxCycleNumber = 0;
        Smin = 0.9;
        Smax = str2double(cycleParameters(2));
        supersaturationLimits = [Smin Smax];
        m21min = initialm21;
        m21max = yieldFactor*m21min;
        cmax = initialConcentration;
        cmin = cmax - particleDensity*shapeFactor*(m21max - m21min);
        Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
        Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

        additiveRange = [0.08 0];
        for k = 1:length(additiveRange)
            additiveConcentration = additiveRange(k);
            solubilityFactor = interp1(additiveFactors(1,:),additiveFactors(2,:),additiveRange(k)); % increase solubility appropriately
            growthFactor = interp1(additiveFactors(1,:),additiveFactors(3,:),additiveRange(k)); % slows down rate of radius growth appropriately

            Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
            Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

% memory efficient high resolution solution:
            [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
            averageRadius = m31./m21;
            averageHeight = m22./m21;
            stdRadius = sqrt(abs(m41./m21 - averageRadius.^2));
            stdHeight = sqrt(abs(m23./m21 - averageHeight.^2));
    
% aspect ratio plots:
            figure(1)
            if k == 1
                tiles1 = tiledlayout(1,2,'Padding','compact','TileSpacing','compact','Units','inches');
                tiles1.InnerPosition(1:2) = [3 3];
                tiles1.InnerPosition(4) = 3.25;
                tiles1.OuterPosition(3) = 7;
                nexttile(1)
                contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.5); % initial PSSD
            end
            nexttile(1)
            hold on
            plot(averageRadius, averageHeight,"LineWidth",1)
            
            % iso-volume plots:
            l1 = 100:0.1:1000;
            avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
            l2 = avgVolumemin./l1.^2;
            plot(l1,l2,'--',"LineWidth",1) % lower bound
            l2 = yieldFactor*avgVolumemin./l1.^2;
            plot(l1,l2,'--',"LineWidth",1) % desired final volume
    
            % ideal aspect ratio:
            plot(L2,L2,':',"LineWidth",1)
            xlim([200 600]);
            ylim([200 800]);
    
            % PSSD plots
            contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.5); % final PSSD
            text(averageRadius(end),averageHeight(end),{['    {\it c}_A=',num2str(additiveConcentration),' [kg kg^{-1}]']},'FontSize',10)
            
            if k == 2
                set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
                colororder(gca,['#000000'; '#d00000';'#03045e'; '#808080'])
                xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',10)
                xticks(linspace(200,600,3))
                ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',10)
                yticks(linspace(200,800,4))
                text(200,600,' {\it V_{min}}','FontSize',10,'Color','#d00000')
                text(320,735,'{\it V_{max}}','FontSize',10,'Color','#03045e')
                text(460,460,'  {\it L_1 = L_2}','FontSize',10,'Color','#808080')
                text(averageRadius(1),averageHeight(1),{'','    Initial','    PSSD'},'FontSize',10)
            end
        end
        

        supersaturationRange = [1.1 1.2];
        for k = 1:length(supersaturationRange)
            Smax = supersaturationRange(k);
            supersaturationLimits = [Smin Smax];

% memory efficient high resolution solution:
            [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
            averageRadius = m31./m21;
            averageHeight = m22./m21;
            stdRadius = sqrt(abs(m41./m21 - averageRadius.^2));
            stdHeight = sqrt(abs(m23./m21 - averageHeight.^2));
    
% aspect ratio plots:
            figure(1)
            if k == 1
                nexttile(2)
                contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.5); % initial PSSD
            end
            nexttile(2)
            hold on
            plot(averageRadius, averageHeight,"LineWidth",1)
            
            % iso-volume plots:
            l1 = 100:0.1:1000;
            avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
            l2 = avgVolumemin./l1.^2;
            plot(l1,l2,'--',"LineWidth",1) % lower bound
            l2 = yieldFactor*avgVolumemin./l1.^2;
            plot(l1,l2,'--',"LineWidth",1) % desired final volume
    
            % ideal aspect ratio:
            plot(L2,L2,':',"LineWidth",1)
            xlim([200 600]);
            ylim([200 800]);
    
            % PSSD plots
            contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.5); % final PSSD
            text(averageRadius(end),averageHeight(end),{['     S=',num2str(Smax)]},'FontSize',10)
            
            if k == 2
                set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left')
                colororder(gca,['#000000'; '#d00000';'#03045e'; '#808080'])
                xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',10)
                xticks(linspace(200,600,3))
                ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',10)
                yticks(linspace(200,800,4))
                text(200,600,' {\it V_{min}}','FontSize',10,'Color','#d00000')
                text(320,735,'{\it V_{max}}','FontSize',10,'Color','#03045e')
                text(460,460,'  {\it L_1 = L_2}','FontSize',10,'Color','#808080')
                text(averageRadius(1),averageHeight(1),{'','    Initial','    PSSD'},'FontSize',10)
                
                fig = gcf;
                fig.WindowState = "maximized";

                fileName = 'growth variation.pdf';
                fpath = strcat(path,'Batch2D_additive\Attainable region growth\', folderName);
                exportgraphics(tiles1, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');
            end
        end

        for k = 1:length(additiveRange)
            additiveConcentration = additiveRange(k);
            solubilityFactor = interp1(additiveFactors(1,:),additiveFactors(2,:),additiveRange(k)); % increase solubility appropriately
            growthFactor = interp1(additiveFactors(1,:),additiveFactors(3,:),additiveRange(k)); % slows down rate of radius growth appropriately

            Smax = supersaturationRange(k);
            supersaturationLimits = [Smin Smax];

            Tmax = (1/0.036)*log(cmax/(3.37*solubilityFactor));
            Tmin = (1/0.036)*log(cmin/(3.37*solubilityFactor));

            % memory efficient high resolution solution:
            [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
            averageRadius = m31./m21;
            averageHeight = m22./m21;
            stdRadius = sqrt(abs(m41./m21 - averageRadius.^2));
            stdHeight = sqrt(abs(m23./m21 - averageHeight.^2));
            
            if k == 1
                figure(2)
                plot(averageRadius, averageHeight,"LineWidth",1)
                hold on
                
                % iso-volume plots:
                l1 = 100:0.1:1000;
                avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
                l2 = avgVolumemin./l1.^2;
                plot(l1,l2,'--',"LineWidth",1) % lower bound
                l2 = yieldFactor*avgVolumemin./l1.^2;
                plot(l1,l2,'--',"LineWidth",1) % desired final volume
        
                % ideal aspect ratio:
                plot(L2,L2,':',"LineWidth",1)
                xlim([200 600]);
                ylim([200 800]);
        
                % PSSD plots
                contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.5); % initial PSSD
                contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.5); % final PSSD
                text(averageRadius(end),averageHeight(end),{'     Lower bound'},'FontSize',10)
                
                set(gca,'FontSize',8,'FontName','times','TitleHorizontalAlignment','left','Units','inches')
                ax = gca;
                ax.InnerPosition(1:2) = [3 3];
                ax.InnerPosition(3:4) = [2.8 3];
                colororder(gca,['#000000'; '#d00000';'#03045e'; '#808080'])
                xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',10)
                xticks(linspace(200,600,3))
                ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',10)
                yticks(linspace(200,800,4))
                text(200,600,' {\it V_{min}}','FontSize',10,'Color','#d00000')
                text(320,735,'{\it V_{max}}','FontSize',10,'Color','#03045e')
                text(350,350,'  {\it L_1 = L_2}','FontSize',10,'Color','#808080')
                text(averageRadius(1),averageHeight(1),{'','    Initial','    PSSD'},'FontSize',10)
                
                fig = gcf;
                fig.WindowState = "maximized";
            end
            if k == 2
                figure(2)
                hold on
                plot(averageRadius, averageHeight,"LineWidth",1)
                contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.5); % final PSSD
                text(averageRadius(end),averageHeight(end),{'     Upper bound'},'FontSize',10)

                fileName = 'growth bounds.pdf';
                fpath = strcat(path,'Batch2D_additive\Attainable region growth\', folderName);
                exportgraphics(gca, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');
            end
        end
        fileName = 'simulation variables.mat';
        fpath = strcat(path,'Batch2D_additive\Attainable region growth\', folderName);
        save(fullfile(fpath, fileName))
        %% Powerpoint figures
    case 9
 % create a folder with the current date and time:
        mkdir(strcat(path,'Batch2D_additive\Powerpoint figures'), folderName)

% prompt user for desired volume and supersaturation limits
        prompt = {'Desired final yield (as a multiple of the original mass):','Number of growth/dissolution cycles:','Minimum supersaturation:','Maximum supersaturation:'};
        dlgtitle = 'Temperature cycle parameters';
        definput = {'3','2','0.9','1.15'};
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

% memory efficient high resolution solution:
        [finalPSD, concentration, ~, ~, ~, ~, m31, m22, m21, ~, ~, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
        averageRadius = m31./m21;
        averageHeight = m22./m21;
        temperatureRamp = [t; temperature]; % to be used by mom
        
% method of moments solution:
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t_mom, y] = ode113(@(t_mom, y)mom2D_additiveTD(t_mom,y,kg11,kg12,kg13,kg21,kg22,kg23,kd11,kd12,kd21,kd22,temperatureRamp,shapeFactor,particleDensity,solubilityFactor,growthFactor), [0 simulationTime], y0,options);
        concentration_mom = y(:,17);
        m00_mom = y(:,1);
        m21_mom = y(:,7);
        m31_mom = y(:,11);
        m22_mom = y(:,9);
        m41_mom = y(:,16);
        m23_mom = y(:,15);
        averageRadius_mom = m31_mom./m21_mom;
        averageHeight_mom = m22_mom./m21_mom;
        stdRadius_mom = sqrt(abs(m41_mom./m21_mom - averageRadius_mom.^2));
        stdHeight_mom = sqrt(abs(m23_mom./m21_mom - averageHeight_mom.^2));

% aspect ratio plots:
        figure(1)
        tiles1 = tiledlayout(1,2,'Padding','compact','TileSpacing','compact','Units','inches');
        tiles1.InnerPosition(1:2) = [2 2];
        tiles1.InnerPosition(4) = 7.2;
        tiles1.OuterPosition(3) = 14;
        nexttile
        plot(averageRadius, averageHeight,"LineWidth",2.5)
        hold on
        plot(averageRadius_mom, averageHeight_mom, '--',"LineWidth",2.5)
        
        % iso-volume plots:
        l1 = 100:0.1:1000;
        avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
        l2 = avgVolumemin./l1.^2;
        plot(l1,l2,'--',"LineWidth",2.5) % lower bound
        l2 = yieldFactor*avgVolumemin./l1.^2;
        plot(l1,l2,'--',"LineWidth",2.5) % desired final volume

        % ideal aspect ratio:
        plot(L2,L2,':',"LineWidth",2.5)
        xlim([200 600]);
        ylim([200 800]);

        % PSSD plots
        contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.5); % initial PSSD
        contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.5); % final PSSD


        set(gca,'FontSize',24,'FontName','times','TitleHorizontalAlignment','left')
        colororder(gca,['#000000'; '#00A300';'#d00000';'#03045e'; '#808080'])
        xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',30)
        xticks(linspace(200,600,3))
        ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',30)
        yticks(linspace(200,800,4))
        % pbaspect([1 1.1 1])
        text(200,600,' {\it V_{min}}','FontSize',24,'Color','#d00000')
        text(320,735,'{\it V_{max}}','FontSize',24,'Color','#03045e')
        text(460,460,'  {\it L_1 = L_2}','FontSize',24,'Color','#808080')
        text(300,525,{'\color{black}FVM','\color[rgb]{0,0.64,0}MOM'},'FontSize',24)
        %legend('High resolution solution','Method of moments solution','Location','southeast')
        text(averageRadius(1),averageHeight(1),{'    Initial','    PSSD'},'FontSize',24)
        text(averageRadius(end),averageHeight(end),{'    Final','    PSSD'},'FontSize',24)
        
% solubility plot:
        nexttile
        colororder(gca,{'#000080','#005000'})
        yyaxis left
        plot(t(temperature<=20),temperature(temperature<=20),'LineWidth',2.5)
        xlim([0 300])
        ylim([-5 35])
        ylabel(['{\it T} [' char(176) 'C]'])
        xlabel('{\it t} [h]')
        yyaxis right
        plot(t,concentration,'LineWidth',2.5)
        ylabel('{\it c} [g kg^{-1}]')
        ylim([3.5 8])
        
        set(gca,'FontSize',24,'FontName','times','TitleHorizontalAlignment','left')

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('aspect ratio + profiles base case.emf');
        fpath = strcat(path,'Batch2D_additive\Powerpoint figures\', folderName);
        exportgraphics(tiles1, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

% Additive concentration plot
        additiveConcentration = 0.08;
        solubilityFactor = interp1(additiveFactors(1,:),additiveFactors(2,:),additiveConcentration); % increase solubility appropriately
        growthFactor = interp1(additiveFactors(1,:),additiveFactors(3,:),additiveConcentration); % slows down rate of radius growth appropriately

        [finalPSD, concentration, G1, G2, supersaturation, m00, m31, m22, m21, m41, m23, t, temperature] = highRes2D_additiveCS(dL1,dL2, L1, L2, simulationTime, kg11, kg12, kg13, kg21, kg22, kg23, kd11, kd12, kd21, kd22, shapeFactor, particleDensity, initialConcentration, initialPSD, yieldFactor, maxCycleNumber, supersaturationLimits,growthFactor,solubilityFactor);
        averageRadius = m31./m21;
        averageHeight = m22./m21;
        stdRadius = sqrt(abs(m41./m21 - averageRadius.^2));
        stdHeight = sqrt(abs(m23./m21 - averageHeight.^2));

% Additive effect
        figure(2)
        plot(averageRadius, averageHeight,"LineWidth",2.5)
        hold on

        % iso-volume plots:
        l1 = 100:0.1:1000;
        avgVolumemin = sum(L1.^4.*L2'.^2.*initialPSD,'all')*dL1*dL2/initialm21; %m42/m21 (volume-weighted average volume)
        l2 = avgVolumemin./l1.^2;
        plot(l1,l2,'--',"LineWidth",2.5) % lower bound
        l2 = yieldFactor*avgVolumemin./l1.^2;
        plot(l1,l2,'--',"LineWidth",2.5) % desired final volume

        % ideal aspect ratio:
        plot(L2,L2,':',"LineWidth",2.5)
        xlim([200 600]);
        ylim([200 800]);

        % PSSD plots:
        contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9],'FaceAlpha',0.6); % initial PSSD
        contourf(L1,L2,finalPSD./max(finalPSD(:)), [0.1 0.5 0.9],'FaceAlpha',0.6); % final PSSD

        set(gca,'FontSize',24,'FontName','times','TitleHorizontalAlignment','left','Units','inches','Position',[2 2 5.5 7.2])
        colororder(['#000000';'#d00000';'#03045e'; '#808080'])
        xlabel({['{\it L}_{1,V} [' char(181) 'm]']},'FontSize',30)
        xticks(linspace(200,600,3))
        ylabel({['{\it L}_{2,V} [' char(181) 'm]']},'FontSize',30)
        yticks(linspace(200,800,4))
        text(200,600,' {\it V_{min}}','FontSize',24,'Color','#d00000')
        text(320,735,'{\it V_{max}}','FontSize',24,'Color','#03045e')
        text(500,500,'  {\it L_1 = L_2}','FontSize',24)
        %legend('High resolution solution','Minimum volume','Desired volume','Ideal aspect ratio','Location','southeast')

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('additive effect (0.08 [kg per kg]).emf');
        fpath = strcat(path,'Batch2D_additive\Powerpoint figures\', folderName);
        exportgraphics(gcf, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');
% PSSD example
        figure(3)
        surf(L1,L2,initialPSD./max(initialPSD(:)),'EdgeColor','none','FaceAlpha',0.8)
        hold on
        contourf(L1,L2,initialPSD./max(initialPSD(:)),[0.1 0.5 0.9]);

        colormap turbo
        set(gca,'FontSize',24,'XTick',[],'YTick',[],'ZTick',[],'fontname','times','Units','inches','Position',[2 2 5.5 7.2])
        xlabel('{\it L}_1','FontSize',30)
        ylabel('{\it L}_2','FontSize',30)
        zlabel('PSSD, {\it f}','FontSize',30)
        xlim([200 300])
        ylim([300 500])

        fig = gcf;
        fig.Units = "normalized";
        fig.OuterPosition = [0 0 1 1];

        fileName = strcat('PSSD surface example.png');
        fpath = strcat(path,'Batch2D_additive\Powerpoint figures\', folderName);
        exportgraphics(gcf, fullfile(fpath, fileName), 'ContentType','image','Resolution',600,'BackgroundColor','none');

        fileName = strcat('PSSD surface example.emf');
        fpath = strcat(path,'Batch2D_additive\Powerpoint figures\', folderName);
        exportgraphics(gcf, fullfile(fpath, fileName), 'ContentType','vector','Resolution',600,'BackgroundColor','none');

% store operating conditions etc in a matrix to be exported as a
% text file
        namesCell = {'Initial Concentration'; 'Additive concentration';'Shape factor'; 'Particle density'; 'Simulation time';'dL1';'dL2';'Yield';'Number of cycles';'Minimum supersaturation';'Maximum supersaturation'};
        unitsCell = {'[g/kg]';'[kg/kg]'; '[-]';['[g/' char(181) 'm' char(179) ']'];'[h]';['[' char(181) 'm]'];['[' char(181) 'm]'];'[-]';'[-]';'[-]';'[-]'};
        parametersCell = cat(1, [namesCell, [essentialParameters(:); cycleParameters], unitsCell]);
        fileName='Iso-volume base case parameters';
        fpath = strcat(path,'Batch2D_additive\Powerpoint figures\', folderName);
        writecell(parametersCell,fullfile(fpath,fileName))

end