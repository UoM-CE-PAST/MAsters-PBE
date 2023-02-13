clear; clc

%% Essential parameters
T=25; %Celsius
ce=3.37*exp(0.0359*T); %g/kg
c0=15; %g/kg
k1 = 10; %um/h
k2 = 1;

%% PSD initialization

%Length range and length step
delL=1; %um
L=1:delL:500; %um

%PSD at time 0
f(:,1)=1e5*normpdf(L,150,20);

%% Finite difference

% Constant growth rate computation
c=c0;
S=c0/ce;
G=5;

% Varying time step
t=0;
i=1;
tmax=30; %h range of t required

while t(i)<tmax
    %Calculate stable time step using available values
    delt=0.9*abs(delL/G);

    %Check if the time step will exceed max value
    if tmax-t(i)<=delt
        t(i+1)=tmax;
        delt=tmax-t(i);
        break
    else
        t(i+1)=t(i)+delt;
    end

    %Calculate the PSD at the new time using backward difference
    j=length(L)-1;
    while j>1
        f(j,i+1)=f(j,i)-G*(delt/delL)*(f(j,i)-f(j-1,i));
        j=j-1;
    end

    % Increase time counter
    i=i+1;
end

%PSD plot

for ii = 1:size(f,2)
    plot(L,f(:,ii))
    pause(1e-4)
end
