function [f, c, G, S, m3, t] = LaxWendroff(delL, L, tmax, k1, k2, kv, T, rhoc, c0, f0)

%Equilibrium conc.
ce=3.37*exp(0.0359*T); % g/kg

%Initial values
f(:,1)=f0;
m3(1)=trapz(L.^3.*f(:,1)');
c(1)=c0;
S(1)=c0/ce;
G(1)=k1*(S(1)-1)^k2;

t=0;
i=1;
while t(i)<tmax
    %Calculate stable time step using available values
    delt=0.9*delL/G(i);

    %Check if the max stable time step will exceed time range
    if tmax-t(i)<=delt
        t(i+1)=tmax;
        delt=tmax-t(i);
    else
        t(i+1)=t(i)+delt;
    end

    %Calculate the PSD at the new time using backwards time difference
    j=length(L)-1;
    while j>2
        f(j,i+1)=f(j,i)-G(i)*0.5*delt/delL*(f(j+1,i)-f(j-1,i))+0.5*(G(i)*delt/delL)^2*(f(j+1,i)-2*f(j,i)+f(j-1,i));
        j=j-1;
    end

    % Calculate essential values for the next time step
    m3(i+1)=trapz(L.^3.*f(:,i+1)');
    c(i+1)=c(i)-rhoc*kv*(m3(i+1)-m3(i));
    S(i+1)=c(i+1)/ce;

    if S(i+1)<=1 % Necessary to make sure it remains a growth problem
        S(i+1)=1;
    end

    G(i+1)=k1*(S(i+1)-1)^k2;
 
    % Increase time counter
    i=i+1;
end
end
