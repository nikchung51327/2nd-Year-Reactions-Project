function [r1, r2] = rates_LHHW(T, alpha, P_CH3OH, P_O2, P_H20, P_HCHO)

    R = 8.314; % J/mol/K
    
    % --- Activation Energy parameters (converted kJ/mol -> J/mol) --- %
    A_kCH3OH = 1.50e7;   Ea_kCH3OH = 86.00e3;
    A_KCH3OH = 2.60e-4;  Ea_KCH3OH = -56.78e3;
    A_KO2    = 1.42e-5;  Ea_KO2    = -60.32e3;
    A_KH2O   = 5.50e-7;  Ea_KH2O   = -86.45e3;
    A_kCO    = 3.50e2;   Ea_kCO    = 46.00e3;
    
    % --- Calculation of rate constants from params --- %
    kCH3OH = A_kCH3OH * exp(-Ea_kCH3OH/(R*T));
    KCH3OH = A_KCH3OH * exp(-Ea_KCH3OH/(R*T));
    KO2    = A_KO2    * exp(-Ea_KO2/(R*T));
    KH2O   = A_KH2O   * exp(-Ea_KH2O/(R*T));
    kCO    = A_kCO    * exp(-Ea_kCO/(R*T));
    
    % --- R1
    nom1 = alpha*kCH3OH*KCH3OH*P_CH3OH*KO2*(sqrt(P_O2));
    den1 = (1+KCH3OH*P_CH3OH + KH2O*P_H20) + (1 + KO2*(sqrt(P_O2)));
    r1 = nom1 / den1;

    % --- R2
    nom2 = alpha*kCO*P_HCHO*KO2*(sqrt(P_O2));
    den2 = 1 + KO2*(sqrt(P_O2));
    r2 = nom2 / den2;

end
