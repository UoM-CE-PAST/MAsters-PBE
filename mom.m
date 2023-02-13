function ddt = mom(y, ce, k1, k2, rhoc, kv)
%% Store variables to be solved

m0 = y(1);
m1 = y(2);
m2 = y(3);
m3 = y(4);
c = y(5);

%% Necessary parameters

S = c/ce;
G = k1*(S-1)^k2;

%% ODE system

dm0dt = 0; % This is added to double check the simulation is correct
dm1dt = G*m0;
dm2dt = 2*G*m1;
dm3dt = 3*G*m2;
dcdt = -rhoc*kv*dm3dt;

ddt = [dm0dt; dm1dt; dm2dt; dm3dt; dcdt];
end