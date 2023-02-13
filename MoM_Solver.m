clear; clc

%% Initial PSD (Gaussian distribution)

L=0:1000;
f0=1e5*normpdf(L,150,20);

%% Numerical integration (to find initial moments)
% step size
h=1;
% Initialization
m0i=trapz(f0);
m1i=trapz(L.*f0);
m2i=trapz(L.^2.*f0);
m3i=trapz(L.^3.*f0);

%% Constant values
T=25; % ambient T
ce = 3.37*exp(0.0359*T); % equilibrium concentration in g/kg
c0 = 30; % g/kg
k1 = 10; % micrometers per s
k2 = 1; % No units
kv = 1; % assuming cuboidal particles
rhoc = 1.46e-12; % g/um3

%% ODE system integration

y0 = [m0i m1i m2i m3i c0];
tspan = [0 100];
[t, y] = ode45(@(t, y)mom(t, y, ce, k1, k2, m0i, rhoc, kv), tspan, y0);

%% Plot of moments

% subplot(3,2,1)
% plot(t,y(:,1))
% 
% subplot(3,2,2)
% plot(t,y(:,2))
% 
% subplot(3,2,3)
% plot(t,y(:,3))
% 
% subplot(3,2,4)
% plot(t,y(:,4))
% 
% subplot(3,2,5)
avgL=y(:,2)./y(:,1);

figure(1)
hold on
subplot(3,1,2)
plot(t,y(:,5),'--','linewidth',1.2), set(gca,'FontSize',18)
legend('Finite Difference','Method of Moments')

hold on
subplot(3,1,3)
plot(t,avgL,'--','linewidth',1.2), set(gca,'FontSize',18)
legend('Finite difference','Method of moments')

