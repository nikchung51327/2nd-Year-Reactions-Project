clc; clear; 
function initial_conditions = get_initial_conditions(m_flowrate, P0, T0)

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

    % --- Total molar flowrate [mol/s] --- %
    Mr_mix = y0_CH3OH*Mr_CH3OH + y0_O2*Mr_O2 + y0_H2O*Mr_H2O + y0_N2*Mr_N2;
    nT0 = m_flowrate / Mr_mix;

    % --- Individual molar flowrates [mol/s] --- %
    n0_CH3OH = y0_CH3OH * nT0;
    n0_O2    = y0_O2    * nT0;
    n0_H2O   = y0_H2O   * nT0;
    n0_HCHO  = y0_HCHO  * nT0;
    n0_CO    = y0_CO    * nT0;
    n0_N2    = y0_N2    * nT0;

    % --- Initial condition vector --- %
    % [CH3OH, O2, HCHO, H2O, CO, N2, extent1, extent2, P, T]
    initial_conditions = [n0_CH3OH, n0_O2, n0_HCHO, n0_H2O, n0_CO, n0_N2, 0, 0, P0, T0];

end

function outlet_conditions = get_outlet_conditions(z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T)

    F_CH3OH_out  = F_CH3OH(end);
    F_O2_out     = F_O2(end);
    F_HCHO_out   = F_HCHO(end);
    F_H2O_out    = F_H2O(end);
    F_CO_out     = F_CO(end);
    F_N2_out     = F_N2(end);
    P_out        = P(end);
    T_out        = T(end);
    extent1_out  = extent1(end);
    extent2_out  = extent2(end);
    outlet_conditions = [F_CH3OH_out, F_O2_out, F_HCHO_out, F_H2O_out, F_CO_out, F_N2_out, extent1_out, extent2_out, P_out, T_out];
end

function [value, isterminal, direction] = T_max_event(z, F, T_max)
    value      = F(10) - T_max;
    isterminal = 1;
    direction  = 1;
end

function [value, isterminal, direction] = P_min_event(z, F, P_min)
    value      = F(9) - P_min;
    isterminal = 1;
    direction  = -1;
end

function [r1, r2] = rates_LHHW(T, alpha, P_CH3OH, P_O2, P_H20, P_HCHO)

    R = 8.314;
    A_kCH3OH = 1.50e7;   Ea_kCH3OH = 86.00e3;
    A_KCH3OH = 2.60e-4;  Ea_KCH3OH = -56.78e3;
    A_KO2    = 1.42e-5;  Ea_KO2    = -60.32e3;
    A_KH2O   = 5.50e-7;  Ea_KH2O   = -86.45e3;
    A_kCO    = 3.50e2;   Ea_kCO    = 46.00e3;

    kCH3OH = A_kCH3OH * exp(-Ea_kCH3OH/(R*T));
    KCH3OH = A_KCH3OH * exp(-Ea_KCH3OH/(R*T));
    KO2    = A_KO2    * exp(-Ea_KO2/(R*T));
    KH2O   = A_KH2O   * exp(-Ea_KH2O/(R*T));
    kCO    = A_kCO    * exp(-Ea_kCO/(R*T));

    P_O2 = sqrt(max(P_O2, 0));
    nom1 = alpha*kCH3OH*KCH3OH*P_CH3OH*KO2*(P_O2);
    den1 = (1+KCH3OH*P_CH3OH + KH2O*P_H20)*(1 + KO2*(P_O2));
    r1 = nom1 / den1;

    nom2 = alpha*kCO*P_HCHO*KO2*(P_O2);
    den2 = 1 + KO2*(P_O2);
    r2 = nom2 / den2;
end

function [r1, r2] = rates_powerlaw(T, P_CH3OH, P_O2, P_HCHO)

    R = 8.314;
    Ea1 = 24.8840e3; Ea2 = 4.4067e4;
    k01 = 1.4614;    k02 = 0.5583;
    k1 = k01 * exp(Ea1/(R*T));
    k2 = k02 * exp(Ea2/(R*T));
    a=0.53; b=0.19; c=0.8635; d=0.2149;
    r1 = k1*(P_CH3OH^a)*(P_O2^b);
    r2 = k2*(P_HCHO^c)*(P_O2^d);
end

function dT_dz = energy_bal(r1, r2, y, n_T, T)

    Hr1 = -159*10^3;
    Hr2 = -200*10^3;
    tube_diameter = 20 * 10^-3;
    A = pi*(tube_diameter/2)^2;

    a = [40.0460; 29.5260; 33.9330; 29.3420; 29.6750; 29.5560];
    b = [-3.8287e-02; -8.8999e-03; -8.4186e-03; -3.5395e-03;  1.8937e-02; -6.5807e-03];
    c = [ 2.4529e-04;  3.8083e-05;  2.9906e-05;  1.0076e-05;  2.8739e-05;  2.0130e-05];
    d = [-2.1679e-07; -3.2629e-08; -1.7825e-08; -4.3116e-09; -2.0092e-08; -1.2227e-08];
    e = [ 5.9909e-11;  8.8607e-12;  3.6934e-12;  2.5935e-13;  3.4333e-12;  2.2617e-12];

    y_Cp = zeros(6,1);
    for i = 1:6
        y_Cp(i) = y(i)*(a(i) + b(i)*T^1 + c(i)*T^2 + d(i)*T^3 + e(i)*T^4);
    end

    sum_y_Cp = sum(y_Cp);
    num = -(r1*Hr1 + r2*Hr2);
    den = n_T*sum_y_Cp;
    dT_dz = num / den;
end

function dni_dz = reactor(t, F, options)

    if F(9) <= 0.01
        fprintf('WARNING: Pressure too low at z=%.2f, P=%.4f\n', t, F(9));
        dni_dz = zeros(10,1);
        return;
    end
    if any(F(1:6) < 0)
        fprintf('WARNING: Negative flows at z=%.2f\n', t);
        dni_dz = zeros(10,1);
        return;
    end

    rho_cat = 1500;
    bed_voidage = 0.4;
    R = 8.314;
    D_p = 7.04e-3;
    alpha = 0.44;

    Mr_CH3OH = 32.042e-3; Mr_O2 = 32.00e-3; Mr_HCHO = 30.026e-3;
    Mr_H2O = 18.00e-3;   Mr_CO = 28e-3;    Mr_N2 = 28.014e-3;
    tube_diameter = 20e-3;
    A = pi*(tube_diameter/2)^2;

    n_CH3OH = F(1); n_O2 = F(2); n_HCHO = F(3); n_H2O = F(4);
    n_CO = F(5);    n_N2 = F(6);
    P = F(9); P_pascal = P * 101325; T = F(10);

    n_T = n_CH3OH + n_O2 + n_HCHO + n_H2O + n_CO + n_N2;
    y_CH3OH = n_CH3OH/n_T; y_O2 = n_O2/n_T; y_H2O = n_H2O/n_T;
    y_N2 = n_N2/n_T;       y_HCHO = n_HCHO/n_T; y_CO = n_CO/n_T;
    y = [y_CH3OH; y_O2; y_H2O; y_N2; y_HCHO; y_CO];

    P_CH3OH = y_CH3OH*P; P_O2 = y_O2*P; P_H20 = y_H2O*P; P_HCHO = y_HCHO*P;

    Mr_total = y_CH3OH*Mr_CH3OH + y_O2*Mr_O2 + y_H2O*Mr_H2O + y_HCHO*Mr_HCHO + y_CO*Mr_CO + y_N2*Mr_N2;
    m_flowrate = Mr_total*n_T;
    G = m_flowrate/A;
    rho_mix = (P_pascal*Mr_total)/(R*T);
    dw_dz = rho_cat*(1-bed_voidage)*A;

    if options(2) == "Yes"
        dP_dz = -(1.75*(G^2)*(1-bed_voidage))/(D_p*rho_mix*(bed_voidage)^3);
        dP_dz = dP_dz / 101325;
    else
        dP_dz = 0;
    end

    if options(1) == "LHHW"
        [r1, r2] = rates_LHHW(T, alpha, P_CH3OH, P_O2, P_H20, P_HCHO);
    else
        [r1, r2] = rates_powerlaw(T, y_CH3OH, y_O2, y_HCHO, y_H2O);
    end

    if options(3) == "Yes"
        dT_dz = energy_bal(r1, r2, y, n_T, T)*dw_dz;
    else
        dT_dz = 0;
    end

    dn_CH3OH_dz  = -r1*dw_dz;
    dn_O2_dz     = -(r1+r2)/2*dw_dz;
    dn_HCHO_dz   = (r1-r2)*dw_dz;
    dn_H2O_dz    = (r1+r2)*dw_dz;
    dn_CO_dz     = r2*dw_dz;
    dn_N2_dz     = 0;
    d_extent1_dz = A*rho_cat*(1-bed_voidage)*r1;
    d_extent2_dz = A*rho_cat*(1-bed_voidage)*r2;

    dni_dz = [dn_CH3OH_dz; dn_O2_dz; dn_HCHO_dz; dn_H2O_dz; dn_CO_dz; dn_N2_dz; d_extent1_dz; d_extent2_dz; dP_dz; dT_dz];
end

function [z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
         solve_reactor(initial_conditions, reactor_length, P0, T0, T_max, T_max_active, P_min, P_min_active, options)

    active_events = {};
    if T_max_active == "Yes"
        active_events{end+1} = @(z,F) T_max_event(z, F, T_max);
    end
    if P_min_active == "Yes"
        active_events{end+1} = @(z,F) P_min_event(z, F, P_min);
    end

    if ~isempty(active_events)
        combined_event = @(z,F) combine_events(z, F, active_events);
        opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'Events', combined_event);
    else
        opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    end

    [z_points, solution, z_stop] = ode15s(@(z, F) reactor(z, F, options), ...
                                           [0 reactor_length], initial_conditions, opts);

    if ~isempty(z_stop)
        if T_max_active == "Yes" && solution(end,10) >= T_max
            fprintf('Reactor stopped at z = %.4f m (T_max = %g K reached)\n', z_stop(end), T_max);
        end
        if P_min_active == "Yes" && solution(end,9) <= P_min
            fprintf('Reactor stopped at z = %.4f m (P_min = %g atm reached)\n', z_stop(end), P_min);
        end
    end

    solution = solution';
    F_CH3OH = solution(1, :);
    F_O2    = solution(2, :);
    F_HCHO  = solution(3, :);
    F_H2O   = solution(4, :);
    F_CO    = solution(5, :);
    F_N2    = solution(6, :);
    extent1 = solution(7, :);
    extent2 = solution(8, :);
    P       = solution(9, :);
    T       = solution(10,:);
end

function [value, isterminal, direction] = combine_events(z, F, event_list)
    value = []; isterminal = []; direction = [];
    for k = 1:length(event_list)
        [v, it, d] = event_list{k}(z, F);
        value      = [value;      v(:)];
        isterminal = [isterminal; it(:)];
        direction  = [direction;  d(:)];
    end
end

function plot_reactor_results(z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T, F_CH3OH_inlet)

    z_points = z_points(:); F_CH3OH = F_CH3OH(:); F_O2 = F_O2(:);
    F_HCHO = F_HCHO(:);     F_H2O = F_H2O(:);     F_CO = F_CO(:);
    F_N2 = F_N2(:);         extent1 = extent1(:);  extent2 = extent2(:);
    P = P(:);                T = T(:);

    figure("Name","Molar Flowrate of Species Across the Reactor");
    plot(z_points, F_CH3OH, 'b-', 'LineWidth', 2); hold on;
    plot(z_points, F_O2,    'r-', 'LineWidth', 2);
    plot(z_points, F_HCHO,  'g-', 'LineWidth', 2);
    plot(z_points, F_H2O,   'c-', 'LineWidth', 2);
    plot(z_points, F_CO,    'm-', 'LineWidth', 2);
    xlabel('Reactor Length (m)'); ylabel('Molar Flow Rate (mol/s)');
    xlim([0 z_points(end)]); legend('CH3OH','O2','HCHO','H2O','CO');
    title("Molar Flowrate of Species Across the Reactor")
    set(gcf,'color','w'); grid on

    figure("Name","Pressure Drop Across Reactor")
    plot(z_points, P, 'LineWidth', 2)
    xlabel('Reactor Length (m)'); ylabel('Pressure (atm)');
    xlim([0 z_points(end)]); title("Pressure Drop Across Reactor")
    set(gcf,'color','w'); grid on

    figure("Name","Extent Of Reaction")
    plot(z_points, extent1, 'b-', 'LineWidth', 2); hold on;
    plot(z_points, extent2, 'r-', 'LineWidth', 2);
    title('Extent Of Reaction'); xlim([0 z_points(end)]);
    legend("Reaction 1","Reaction 2"); set(gcf,'color','w'); grid on

    figure("Name","Temperature")
    plot(z_points, T, 'LineWidth', 2)
    title('Temperature Profile'); xlim([0 z_points(end)]);
    set(gcf,'color','w'); grid on

    r1_profile = zeros(length(z_points), 1);
    r2_profile = zeros(length(z_points), 1);
    for i = 1:length(z_points)
        n_T_i     = F_CH3OH(i)+F_O2(i)+F_HCHO(i)+F_H2O(i)+F_CO(i)+F_N2(i);
        y_CH3OH_i = F_CH3OH(i)/n_T_i; y_O2_i = F_O2(i)/n_T_i;
        y_H2O_i   = F_H2O(i)/n_T_i;  y_HCHO_i = F_HCHO(i)/n_T_i;
        P_CH3OH_i = y_CH3OH_i*P(i); P_O2_i = y_O2_i*P(i);
        P_H2O_i   = y_H2O_i*P(i);   P_HCHO_i = y_HCHO_i*P(i);
        [r1_profile(i), r2_profile(i)] = rates_LHHW(T(i), 0.44, P_CH3OH_i, P_O2_i, P_H2O_i, P_HCHO_i);
    end

    figure("Name","Reaction Rates Across Reactor");
    plot(z_points, r1_profile, 'b-', 'LineWidth', 2); hold on;
    plot(z_points, r2_profile, 'r-', 'LineWidth', 2);
    xlabel('Reactor Length (m)'); ylabel('Rate (mol/kg_{cat}/s)');
    legend('r_1 (CH_3OH \rightarrow HCHO)','r_2 (HCHO \rightarrow CO)');
    xlim([0 z_points(end)]); title('Reaction Rates Across the Reactor'); grid on

    X_CH3OH = (F_CH3OH_inlet - F_CH3OH) / F_CH3OH_inlet;
    figure("Name","Conversion of CH3OH Across Reactor")
    plot(z_points, X_CH3OH, 'b-', 'LineWidth', 2)
    xlabel('Reactor Length (m)'); ylabel('Conversion X_{CH_3OH}')
    title('Conversion of CH_3OH Across the Reactor')
    xlim([0 z_points(end)]); ylim([0 1]); set(gcf,'color','w'); grid on
end

% --- Interstage cooling: n_stages defined, each stage runs to T_max boundary --- %
function [z_combined, F_CH3OH_combined, F_O2_combined, F_HCHO_combined, ...
          F_H2O_combined, F_CO_combined, F_N2_combined, ...
          extent1_combined, extent2_combined, P_combined, T_combined, F_CH3OH_inlet] = ...
          reactor_interstage_cooling(m_flowrate, n_stages, T_subcool, P0, T0, T_max, T_max_active, P_min, P_min_active, options)

    % Each stage runs until T_max (or P_min) is hit - no fixed length needed
    max_stage_length = 50;  % [m] safety cap - should never be reached if T_max_active = "Yes"

    IC = get_initial_conditions(m_flowrate, P0, T0);
    F_CH3OH_inlet = IC(1);

    z_combined        = [];
    F_CH3OH_combined  = [];
    F_O2_combined     = [];
    F_HCHO_combined   = [];
    F_H2O_combined    = [];
    F_CO_combined     = [];
    F_N2_combined     = [];
    extent1_combined  = [];
    extent2_combined  = [];
    P_combined        = [];
    T_combined        = [];

    z_offset = 0;

    for s = 1:n_stages

        fprintf('--- Stage %d | T_in = %.1f K ---\n', s, IC(10));

        [z_s, F_CH3OH_s, F_O2_s, F_HCHO_s, F_H2O_s, F_CO_s, F_N2_s, ...
         extent1_s, extent2_s, P_s, T_s, z_stop_s] = ...
            solve_reactor(IC, max_stage_length, P0, T0, T_max, T_max_active, P_min, P_min_active, options);

        fprintf('    Stage %d ended at z = %.4f m | T_out = %.1f K | X = %.3f\n', ...
            s, z_s(end), T_s(end), (F_CH3OH_inlet - F_CH3OH_s(end)) / F_CH3OH_inlet);

        z_combined        = [z_combined;        z_s(:) + z_offset];
        F_CH3OH_combined  = [F_CH3OH_combined,  F_CH3OH_s];
        F_O2_combined     = [F_O2_combined,     F_O2_s];
        F_HCHO_combined   = [F_HCHO_combined,   F_HCHO_s];
        F_H2O_combined    = [F_H2O_combined,    F_H2O_s];
        F_CO_combined     = [F_CO_combined,      F_CO_s];
        F_N2_combined     = [F_N2_combined,     F_N2_s];
        extent1_combined  = [extent1_combined,  extent1_s];
        extent2_combined  = [extent2_combined,  extent2_s];
        P_combined        = [P_combined,        P_s];
        T_combined        = [T_combined,        T_s];

        z_offset = z_offset + z_s(end);

        if s < n_stages
            IC = get_outlet_conditions(z_s, F_CH3OH_s, F_O2_s, F_HCHO_s, F_H2O_s, ...
                                       F_CO_s, F_N2_s, extent1_s, extent2_s, P_s, T_s);
            IC(10) = T_subcool;  % reset temperature to subcooled value
            fprintf('    Cooled to %.1f K for stage %d\n', T_subcool, s+1);
        end
    end
end

%%--------------------------------------------------------- 
%---------------------------------------------------------- 
%------------------ONLY INPUT HERE-------------------------
%----------------------------------------------------------
m_flowrate      = 0.001;
P0              = 1.6;    % Initial Pressure [atm]
T0              = 460;    % Initial Temperature [K]
Kinetic_type    = "LHHW";
Pressure_drop   = "Yes";
Temp_dependence = "Yes";
T_max           = 700;    % Maximum temperature [K]
T_max_active    = "Yes";
P_min           = 1.1;    % Minimum pressure [atm]
P_min_active    = "Yes";
T_subcool       = 670;    % Interstage cooling temperature [K]
n_stages        = 10;      % Number of stages
%----------------------------------------------------------
%----------------------------------------------------------
%----------------------------------------------------------
option = [Kinetic_type Pressure_drop Temp_dependence];

[z, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T, F_CH3OH_inlet] = ...
    reactor_interstage_cooling(m_flowrate, n_stages, T_subcool, P0, T0, ...
                               T_max, T_max_active, P_min, P_min_active, option);

plot_reactor_results(z, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, ...
                     extent1, extent2, P, T, F_CH3OH_inlet)

F_HCHO_outlet = F_HCHO(end);

no_of_tubes = 38.6/F_HCHO_outlet;
fprintf("Number of tubes required: %.0f", no_of_tubes)

% %% colormap for how m_flow and length affects conversion
% 
% T_max_active = "No";
% L_range = 0.5:0.001:2;
% m_flowrange = 0.0005:0.00001:0.001;
% 
% HCHO_conversion = zeros(length(m_flowrange), length(L_range));
% z_stop_heatmap = zeros(length(m_flowrange), 1);  % z_stop using same m_flowrange
% 
% total_iters = length(m_flowrange) * length(L_range);
% current_iter = 0;
% 
% wb = waitbar(0, 'Running simulations...', 'Name', 'Heatmap Progress');
% 
% for i = 1:length(m_flowrange)
%     m_flow = m_flowrange(i);
% 
%     % Get z_stop for this mass flow at full reactor length
%     [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, z_stop_i] = ...
%         solve_reactor(m_flow, reactor_length, P0, T0, T_max, "Yes", option);
%     if isempty(z_stop_i)
%         z_stop_heatmap(i) = reactor_length;  % never hit T_max
%     else
%         z_stop_heatmap(i) = z_stop_i(end);
%     end
% 
%     for j = 1:length(L_range)
%         L = L_range(j); 
%         [z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
%             solve_reactor(m_flow, L, P0, T0, T_max, T_max_active, option);
%         HCHO_conversion(i,j) = F_HCHO(end) / F_CH3OH(1);
% 
%         current_iter = current_iter + 1;
%         waitbar(current_iter / total_iters, wb, ...
%             sprintf('Running: m\\_flow %d/%d, L %d/%d  (%.1f%%)', ...
%             i, length(m_flowrange), j, length(L_range), 100*current_iter/total_iters));
%     end
% end
% 
% close(wb);
% 
% %%
% figure
% imagesc(L_range, m_flowrange, HCHO_conversion)
% set(gca, 'YDir', 'normal')
% xlabel('Reactor Length (m)')
% ylabel('Mass Flow per Tube (kg/s)')
% title('HCHO Yield Heatmap with T_{max} Boundary')
% clim([0 1])
% colorbar
% colormap("turbo")
% hold on
% plot(z_stop_heatmap, m_flowrange, 'w--', 'LineWidth', 2)
% legend('T_{max} boundary', 'Location', 'northwest')
% set(gcf,'color','w')