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
f0 = 1e5*normpdf(L,150,20); % Gaussian

% f0 = zeros(1,length(L)); % Pulse
% f0(100:200)=1e3;

%% Base Case

[f, c, G, S, m3, t] = LaxWendroff(delL, L, tmax, k1, k2, kv, T, rhoc, c0, f0);

%Moving plot
% figure(1)
% subplot(2,1,1)
% for ii = 1:size(f,2)
%     plot(L,f(:,ii), 'linewidth',1.2)
%     pause(1e-5)
% end
% xlabel('Length (um)'), ylabel('f (/um kg)')

%Static plot
subplot(2,1,1)
plot(L,f(:,1),'linewidth',1.2), hold on, plot(L, f(:,end),'linewidth',1.2)
xlabel('Length (um)'), ylabel('f (/um kg)')
legend('Inital PSD','Final PSD')

subplot(2,1,2)
plot(t,c,'linewidth',1.2), xlabel('time (h)'),ylabel('concentration (g/kg)')