function main_rxny2_t3()

    %% --- User Inputs --- %%
    % -- Catalyst info -- %
    % catalyst density [kg/m^3]
    rho_cat = 1500;

    % bed voidage fraction
    bed_voidage = 0.4;

    % tube diameter [mm]; 
    tube_diameter = 20;
    % ------------------- %

    % -- Operating conditions -- %
    % Total operating pressure (assume constant for task3) [atm]
    Ptot = 1.6;
    
    % Temperature [K]
    T = 460;

    % mol fractions; molar mass [kg/mol]
    y0_CH3OH = 0.11; Mr_CH3OH = 0.03204;
    y0_O2 = 0.06; Mr_O2 = 0.032;
    y0_H2O = 0.02; Mr_H2O = 0.01802;
    y0_HCHO = 0; 
    y0_CO = 0;
    y0_N2 = 0.81; Mr_N2 = 0.028;

    % extents of reaction
    extent1 = 0;
    extent2 = 0;

    % finding total flowrate [mol/s] given inlet mass flowrate [kg/s].
    mT0 = 1;
    nT0 = (mT0)/(y0_CH3OH*Mr_CH3OH+y0_O2*Mr_O2+y0_H2O*Mr_H2O+y0_N2*Mr_N2);

    % --------------------- %

    % -- For calculation --%
    
    % LHHW constant alpha
    alpha = 0.44;

    % integrating limits of weight
    Wspan = [0 4];

    % type of kinetics [LHHW, POWER]
    kinetics = 'LHHW';
    
    % Initial Conditions

    F0_CH3OH = y0_CH3OH*nT0; 
    F0_O2 = y0_O2*nT0;
    F0_H2O = y0_H2O*nT0;
    F0_HCHO = y0_HCHO*nT0; 
    F0_CO = y0_CO*nT0;
    F0_N2 = y0_N2*nT0;

    % ------------------- %

    %%% ------------------------------------- %%%
    
    F0_ls = [F0_CH3OH; F0_O2; F0_H2O; F0_HCHO; F0_CO; F0_N2; extent1; extent2];

    
    
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [W, F] = ode45(@(W,f) odes_rxny2_t3(W, f, T, alpha, kinetics, Ptot), Wspan, F0_ls, opts);

    F_CH3OH = F(:,1);
    F_O2 = F(:,2);
    F_H2O = F(:,3);
    F_HCHO = F(:,4);
    F_CO = F(:,5);
    F_N2 = F(:,6);
    extent1 = F(:,7);
    extent2 = F(:,8);
    
    % Tube area [m^2]
    area_tube = pi*((tube_diameter/1000)^2)/4;

    Z = W./(rho_cat*(1-bed_voidage)*area_tube);

    figure; plot(Z,F_CH3OH,Z,F_O2,Z,F_HCHO,Z,F_H2O,Z,F_CO,Z,F_N2,'LineWidth',1.5);
    xlabel('Length of reactor (m^3)'); ylabel('Molar flow F_i (mol/s)');
    legend('CH_3OH','O_2','HCHO','H_2O','CO', 'N2'); grid on; title('Flowrates vs Length');

    figure; plot(W,F_CH3OH,W,F_O2,W,F_HCHO,W,F_H2O,W,F_CO,W,F_N2,'LineWidth',1.5);
    xlabel('Catalyst weight W (kg_{cat})'); ylabel('Molar flow F_i (mol/s)');
    legend('CH_3OH','O_2','HCHO','H_2O','CO', 'N2'); grid on; title('Flowrates vs Weight');

    figure; plot(extent1,Z,'LineWidth',1.5); hold on; plot(extent2,Z,'LineWidth',1.5);
    xlabel('Length of reactor (m^3)'); ylabel('Extent \xi (mol/s)');
    legend('\xi_1 (R1)','\xi_2 (R2)'); grid on; title('Extents vs Length');

    figure; plot(extent1,W,'LineWidth',1.5); hold on; plot(extent2,W,'LineWidth',1.5);
    xlabel('Catalyst weight W (kg_{cat})'); ylabel('Extent \xi (mol/s)');
    legend('\xi_1 (R1)','\xi_2 (R2)'); grid on; title('Extents vs Weight');

end







