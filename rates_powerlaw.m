function [r1, r2] = rates_powerlaw(T, P_CH3OH, P_O2, P_HCHO)
    
    % --- Inputs:
    %
    % T ~ Temperature [K]
    % 
    % P_i ~ Partial pressures of each component [atm]
    %
    % --- Outputs:
    %
    % r1 ~ rate of partial oxidation of ethanol [mol/(kg*s)]
    % r2 ~ rate of partial oxidation of carbon monoxide [mol/(kg*s)]


    R = 8.314; % [J/(mol*K)]
    
    % --- Activation Energy parameters [J/mol] --- %
    
    Ea1 = 24.8840e3; Ea2 = 4.4067e4;
    
    % --- k0 parameters --- %

    k01 = 1.4614; k02 = 0.5583;
    
    % --- Calculation of k1 & k2 --- %
    k1 = k01 * exp(Ea1/(R*T));
    k2 = k02 * exp(Ea2/(R*T));

    % --- Powers for each species --- %

    % r1, orders
    % wrt CH3OH
    a=0.53;

    % wrt O2 
    b=0.19;

    %r2, orders averaged
    % wrt HCHO 
    c=0.8635; 

    % wrt O2
    d=0.2149;

    %P_CH3OH = P_CH3OH*10^5;
    %P_O2 = P_O2*10^5;
    %P_HCHO = P_HCHO*10^5;

    r1 = k1*(P_CH3OH^a)*(P_O2^b);
    r2 = k2*(P_HCHO^c)*(P_O2^d);

end