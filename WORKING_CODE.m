function [r1, r2] = rates_LHHW(T, alpha, P_CH3OH, P_O2, P_H20, P_HCHO)
    
    % --- Inputs:
    %
    % T ~ Temperature [K]
    % alpha ~ Fraction of the mass of catalyst containing the catalytic compounds
    % P_i ~ Partial pressures of each component [atm]
    %
    % --- Outputs:
    %
    % r1 ~ rate of partial oxidation of ethanol [mol/(kg*s)]
    % r2 ~ rate of partial oxidation of carbon monoxide [mol/(kg*s)]

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
    
    P_O2 = sqrt(max(P_O2, 0));
    % --- r1
    nom1 = alpha*kCH3OH*KCH3OH*P_CH3OH*KO2*(P_O2);
    den1 = (1+KCH3OH*P_CH3OH + KH2O*P_H20)*(1 + KO2*(P_O2));
    r1 = nom1 / den1;

    % --- r2
    nom2 = alpha*kCO*P_HCHO*KO2*(P_O2);
    den2 = 1 + KO2*(P_O2);
    r2 = nom2 / den2;
end

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

function dT_dz = energy_bal(r1, r2, y, n_T, T)

    % Enthalpy of reaction at 298K kJ/mol
    Hr1 = -159*10^3;
    Hr2 = -200*10^3;

    % tube diameter [m]; 
    tube_diameter = 20 * 10^-3;

    A = pi*(tube_diameter/2)^2;
    
    species = {'CH3OH'; 'O2'; 'H2O'; 'N2'; 'HCHO'; 'CO'};
    a = [40.0460; 29.5260; 33.9330; 29.3420; 29.6750; 29.5560];
    b = [-3.8287e-02; -8.8999e-03; -8.4186e-03; -3.5395e-03;  1.8937e-02; -6.5807e-03];
    c = [ 2.4529e-04;  3.8083e-05;  2.9906e-05;  1.0076e-05;  2.8739e-05;  2.0130e-05];
    d = [-2.1679e-07; -3.2629e-08; -1.7825e-08; -4.3116e-09; -2.0092e-08; -1.2227e-08];
    e = [ 5.9909e-11;  8.8607e-12;  3.6934e-12;  2.5935e-13;  3.4333e-12;  2.2617e-12];
    
    y_Cp = zeros(6,1); % y*Cp
    
    for i = 1:6
    y_Cp(i) = y(i)*(a(i)* + b(i)*T^1 + c(i)*T^2 + d(i)*T^3 + e(i)*T^4);
    end
    
    sum_y_Cp = sum(y_Cp); 
    num = - (r1*Hr1 + r2*Hr2)*A;
    den = n_T*sum_y_Cp;
    
    dT_dz = num / den;
end

function dni_dz = reactor(t,F,options)


if F(9) <= 0.01  % If pressure too low
    fprintf('WARNING: Pressure too low at z=%.2f, P=%.4f\n', t, F(9));
    dni_dz = zeros(10,1);
    return;
end

% Check for negative molar flows
if any(F(1:6) < 0)
    fprintf('WARNING: Negative flows at z=%.2f\n', t);
    dni_dz = zeros(10,1);
    return;
end


%parameters
rho_cat = 1500; %kg m^-3
bed_voidage = 0.4;
%
R = 8.314;
D_p = 7.04 *10^-3;
%alpha
alpha = 0.44;
% Molar masses of each component kg.mol-1
Mr_CH3OH = 32.042e-3;
Mr_O2 = 32.00e-3;
Mr_HCHO = 30.026e-3;
Mr_H2O = 18.00e-3;
Mr_CO = 28e-3;
Mr_N2 = 28.014e-3;

% tube diameter [m]; 
tube_diameter = 20 * 10^-3;

A = pi*(tube_diameter/2)^2;

%defining each variable 
n_CH3OH = F(1);
n_O2 = F(2);
n_HCHO = F(3);
n_H2O = F(4);
n_CO = F(5);
n_N2 = F(6);
extent_1 = F(7);
extent_2 = F(8);
P = F(9);
P_pascal = P * 101325;% Convert pressure from atm to pascal
T = F(10);

%total molar flowrate
n_T = n_CH3OH + n_O2 + n_HCHO + n_H2O + n_CO + n_N2;

%mol fraction
y_CH3OH = n_CH3OH / n_T;
y_O2 = n_O2 / n_T;
y_H2O = n_H2O / n_T;
y_N2 = n_N2/ n_T;
y_HCHO = n_HCHO / n_T;
y_CO = n_CO/n_T;


y = [y_CH3OH; y_O2; y_H2O; y_N2; y_HCHO; y_CO];
%partial pressures 
P_CH3OH = y_CH3OH * P; 
P_O2 = y_O2 * P; 
P_H20 = y_H2O * P; 
P_HCHO = y_HCHO * P;
P_CO = y_CO * P;
P_N2 = y_N2 * P;

Mr_total = y_CH3OH*Mr_CH3OH + y_O2*Mr_O2 + y_H2O*Mr_H2O + y_HCHO*Mr_HCHO + y_CO*Mr_CO + y_N2*Mr_N2;

m_flowrate = Mr_total*n_T;
%G (superficial mass velocity)
G = m_flowrate/A;  %superficial mass velocity of the fluid
rho_mix = (P_pascal*Mr_total)/(R*T);


dw_dz = rho_cat*(1-bed_voidage)*A;

%Pressure Drop - Yes or No
if options(2) == "Yes"
    %Ergun Equation 
    dP_dz = -(1.75*(G^2)*(1-bed_voidage))/(D_p*rho_mix*(bed_voidage)^3);
    dP_dz = dP_dz / 101325;  % convert Pa/m to atm/m
end 
if options(2) == "No"
    dP_dz = 0;
end 


%rates 
if options(1) == "LHHW"
    [r1, r2] = rates_LHHW(T, alpha, P_CH3OH, P_O2, P_H20, P_HCHO);
end 
if options(1) == "Simple Power"
    [r1, r2] = rates_powerlaw(T, y_CH3OH, y_O2, y_HCHO, y_H2O);
end 

%Isothermal - Yes or No
if options(3) == "Yes"
    dT_dz = energy_bal(r1, r2, y, n_T, T);
end 
if options(3) == "No"
    dT_dz = 0;
end 

dn_CH3OH_dz = -r1*dw_dz;
dn_O2_dz = -(r1+r2)/2*dw_dz;
dn_HCHO_dz = (r1-r2)*dw_dz;
dn_H2O_dz = (r1+r2)*dw_dz; 
dn_CO_dz = r2*dw_dz; 
dn_N2_dz = 0; % Assuming no change in N2 concentration
d_extent1_dz = A*r1;
d_extent2_dz = A*r2;


dni_dz = [dn_CH3OH_dz; dn_O2_dz; dn_HCHO_dz; dn_H2O_dz; dn_CO_dz; dn_N2_dz; d_extent1_dz; d_extent2_dz; dP_dz; dT_dz];


end 
function [z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T] = solve_reactor(m_flowrate,reactor_length,P0,T0,options)


    % --- Inputs:
    % m_flowrate [kg/s]
    % option array
    % Kinetic_type = "LHHW" or "Simple Power";
    % Pressure_drop = "Yes" or "No";
    % Temp_dependence = "Yes" or "No";


    % --- Molar masses [kg/mol] --- %
    Mr_CH3OH = 32.042e-3;
    Mr_O2    = 32.00e-3;
    Mr_H2O   = 18.00e-3;
    Mr_N2    = 28.014e-3;

    % --- Inlet mol fractions --- %
    y0_CH3OH = 0.11;
    y0_O2    = 0.06;
    y0_H2O   = 0.02;
    y0_HCHO  = 0;
    y0_CO    = 0;
    y0_N2    = 0.81;

    % --- Initial molar flow rates [mol/s] --- %
    nT0 = m_flowrate / (y0_CH3OH*Mr_CH3OH + y0_O2*Mr_O2 + y0_H2O*Mr_H2O + y0_N2*Mr_N2);

    n0_CH3OH = y0_CH3OH * nT0;
    n0_O2    = y0_O2    * nT0;
    n0_H2O   = y0_H2O   * nT0;
    n0_HCHO  = y0_HCHO  * nT0;
    n0_CO    = y0_CO    * nT0;
    n0_N2    = y0_N2    * nT0;
    
    %Initial Condition Array
    %1 - CH3OH flowrate
    %2 - O2 flowrate
    %3 - HCHO flowrate
    %4 - H2O flowrate
    %5 - CO flowrate
    %6 - N2 flowrate
    %7 - Extent of reaction 1
    %8 - Extent of reaction 2
    %9 - Pressure 
    %10 - Temperature 
    initial_conditions = [n0_CH3OH, n0_O2, n0_HCHO, n0_H2O, n0_CO, n0_N2, 0, 0, P0, T0];

    % --- Solve ODE --- %
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    [z_points, solution] = ode15s(@(z, F) reactor(z, F, options), [0 reactor_length], initial_conditions, opts);
    solution = solution';

    %output
    F_CH3OH = solution(1, :);
    F_O2    = solution(2, :);
    F_HCHO  = solution(3, :);
    F_H2O   = solution(4, :);
    F_CO    = solution(5, :);
    F_N2    = solution(6, :);
    extent1 = solution(7, :);
    extent2 = solution(8, :);
    P = solution(9,:);
    T = solution(10,:);

end


%---------------------------------------------------------- 
%------------------ONLY INPUT HERE-------------------------
%----------------------------------------------------------
m_flowrate = 0.001;
reactor_length = 2;
P0 = 1.6; %Initial Pressure [bar]
T0 = 523; %Initial Temperature [K]
Kinetic_type = "LHHW";
Pressure_drop = "No";
Temp_dependence = "Yes";
%----------------------------------------------------------
%----------------------------------------------------------
%----------------------------------------------------------

option = [Kinetic_type Pressure_drop Temp_dependence];

[z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T]= solve_reactor(m_flowrate,reactor_length,P0,T0,option);

figure("Name","Molar Flowrate of Species Across the Reactor");
plot(z_points, F_CH3OH, 'b-', 'LineWidth', 2); hold on;
plot(z_points, F_O2, 'r-', 'LineWidth', 2);
plot(z_points, F_HCHO, 'g-', 'LineWidth', 2);
plot(z_points, F_H2O, 'c-', 'LineWidth', 2);
plot(z_points, F_CO, 'm-', 'LineWidth', 2);
plot(z_points, F_N2, 'm-', 'LineWidth', 2);
xlabel('Reactor Length (m)');
ylabel('Molar Flow Rate (mol/s)');
legend('CH3OH', 'O2', 'HCHO', 'H2O', 'CO', 'N2');
title("Molar Flowrate of Species Across the Reactor")

figure("Name","Pressure Drop Across Reactor")
plot(z_points, P,'LineWidth',2)
xlabel('Reactor Length (m)');
ylabel('Pressure (atm)');
title("Pressure Drop Across Reactor")

figure("Name","Extent Of reaction")
plot(z_points, extent1, 'b-', 'LineWidth', 2); hold on;
plot(z_points, extent2, 'r-', 'LineWidth', 2);
title('Extent Of reaction');
legend("Reaction 1","Reaction 2")

figure("Name", "Temperature")
plot(z_points, T, 'LineWidth', 2)
title('Temperature Profile')
legend('T')