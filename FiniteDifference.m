clear; clc

%% Essential parameters
T=25; % C
c0=30; % g/kg
k1 = 10; % um/h
k2 = 1;
kv = 1; % assuming cuboidal particles
rhoc = 1.46e-12; % g/um3
tmax=100; %h range of t required

%Length range and length step
delL=1; %um
L=0:delL:1000; %um

%Initial PSD

% Gaussian
f0 = 1e5*normpdf(L,150,20); 

% Pulse
%f0 = zeros(1,length(L));
% f0(100:200)=1e3; 

%% Base Case

[f, c, G, S, m3, t] = upwind(delL, L, tmax, k1, k2, kv, T, rhoc, c0, f0);

%Moving plot
% figure(1)
% subplot(2,1,1)
% for ii = 1:size(f,2)
%     plot(L,f(:,ii), 'linewidth',1.2)
%     pause(1e-5)
% end
% xlabel('Length (um)'), ylabel('f (/um kg)')

%Static plot
hold on
subplot(3,1,1)
plot(L,f(:,1), 'linewidth',1.2), xlabel('Length [μm]'), ylabel('f [μm^{-1} kg^{-1}]'), hold on, plot(L,f(:,end), 'linewidth',1.2), xlabel('Length [μm]'), ylabel('f [μm^{-1} kg^{-1}]')
legend('Inital PSD','Final PSD')
set(gca,'FontSize',18)
title('(a)','FontSize',24)

hold on
subplot(3,1,2)
plot(t,c,'linewidth',1.2), xlabel('time [h]'),ylabel('concentration [g kg^{-1}]')
set(gca,'FontSize',18)
title('(b)','FontSize',24)

m0 = zeros(length(L),length(t));
m1 = zeros(length(L),length(t));
for i=1:length(t)
m0(1,i) = trapz(f(:,i));
m1(1,i) = trapz(L'.*f(:,i));
end
avgL=m1./m0;

hold on
subplot(3,1,3)
plot(t,avgL(1,:),'LineWidth',1.2), xlabel('time [h]'), ylabel('Average length [μm]')
set(gca,'FontSize',18)
title('(c)','FontSize',24)

%% Temperature variation
% 
% T = [5 15 25 35];
% 
% figure(2)
% for k=1:length(T)
%     clear c t f G S m3
%     [f, c, G, S, m3, t] = upwind(delL, L, tmax, k1, k2, kv, T(k), rhoc, c0, f0);
% 
%     figure(2)
%     subplot(2,1,1)
%     if k==1, plot(L,f(:,1),'--', 'linewidth',1.2), hold on, end
%     plot(L,f(:,end), 'linewidth',1.2), xlabel('Length [μm]'), ylabel('f [μm^{-1} kg^{-1}]'), hold on, set(gca,'FontSize',18)
%     if k==length(T), legend('Initial PSD','T=5 C','T=15 C','T= 25 C','T= 35 C'); end
%     title('(a)','FontSize',24)
% 
%     subplot(2,1,2)
%     plot(t,c,'linewidth',1.2), xlabel('time [h]'),ylabel('concentration [g kg^{-1}]'), hold on, set(gca,'FontSize',18)
%     if k==length(T), legend('T=5 C','T=15 C','T= 25 C','T= 35 C'); end
%     title('(b)','FontSize',24)
% end
% %% Initial concentration variation
% T = 25; % C reset
% c0 = [10 20 30 40 50];
% 
% for k=1:length(c0)
%     clear c t f G S m3
%     [f, c, G, S, m3, t] = upwind(delL, L, tmax, k1, k2, kv, T, rhoc, c0(k), f0);
% 
%     figure(3)
%     subplot(2,1,1)
%     if k==1, plot(L,f(:,1),'--', 'linewidth',1.2), hold on, end
%     plot(L,f(:,end), 'linewidth',1.2), xlabel('Length [μm]'), ylabel('f [μm^{-1} kg^{-1}]'), hold on, set(gca,'FontSize',18)
%     if k==length(c0), legend('Initial PSD','c_{0} = 10 g/kg', 'c_{0} = 20 g/kg', 'c_{0} = 30 g/kg', 'c_{0} = 40 g/kg', 'c_{0} = 50 g/kg'); end
%     title('(a)','FontSize',24)
% 
%     subplot(2,1,2)
%     plot(t,c,'linewidth',1.2), xlabel('time [h]'),ylabel('concentration [g kg^{-1}]'), hold on, set(gca,'FontSize',18)
%     if k==length(c0), legend('c_{0} = 10 g/kg', 'c_{0} = 20 g/kg', 'c_{0} = 30 g/kg', 'c_{0} = 40 g/kg', 'c_{0} = 50 g/kg'); end
%     title('(b)','FontSize',24)
% end
% 
% %% k1 variation
% c0 = 30; %g/kg
% k1 = [5 10 15 20];
% tmax =100;
% for k=1:length(k1)
%     clear c t f G S m3
%     [f, c, G, S, m3, t] = upwind(delL, L, tmax, k1(k), k2, kv, T, rhoc, c0, f0);
% 
%     figure(4)
%     subplot(2,1,1)
%     if k==1, plot(L,f(:,1),'--', 'linewidth',1.2), hold on, end
%     plot(L,f(:,end), 'linewidth',1.2), xlabel('Length [μm]'), ylabel('f [μm^{-1} kg^{-1}]'), hold on, set(gca,'FontSize',18)
%     if k==length(k1), legend('Initial PSD','k_{1} = 5', 'k_{1} = 10', 'k_{1} = 15', 'k_{1} = 20'); end
%     title('(a)','FontSize',24)
% 
%     subplot(2,1,2)
%     plot(t,c,'linewidth',1.2), xlabel('time [h]'),ylabel('concentration [g kg^{-1}]'), hold on, set(gca,'FontSize',18)
%     if k==length(k1), legend('k_{1} = 5', 'k_{1} = 10', 'k_{1} = 15', 'k_{1} = 20'); end
%     title('(b)','FontSize',24)
% end
% 
% %% k2 variation
% k1 = 10;
% k2 = [0.5 1 2];
% 
% tmax = 200;
% for k=1:length(k2)
%     clear c t f G S m3
%     [f, c, G, S, m3, t] = upwind(delL, L, tmax, k1, k2(k), kv, T, rhoc, c0, f0);
% 
%     figure(5)
%     subplot(2,1,1)
%     if k==1, plot(L,f(:,1),'--', 'linewidth',1.2), hold on, end
%     plot(L,f(:,end), 'linewidth',1.2), xlabel('Length [μm]'), ylabel('f [μm^{-1} kg^{-1}]'), hold on, set(gca,'FontSize',18)
%     if k==length(k2), legend('Initial PSD','k_{2} = 0.5', 'k_{2} = 1', 'k_{2} = 2'); end
%     title('(a)','FontSize',24)
% 
%     subplot(2,1,2)
%     plot(t,c,'linewidth',1.2), xlabel('time [h]'),ylabel('concentration [g kg^{-1}]'), hold on, set(gca,'FontSize',18)
%     if k==length(k2), legend('k_{2} = 0.5', 'k_{2} = 1', 'k_{2} = 2'); end
%     title('(b)','FontSize',24)
% end