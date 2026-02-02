function P_i = partial_pressures(F_i, Ptot)

    
    % --- Inputs:
    %
    % F_i ~ Flowrate of each component [mol/s]
    % P_tot ~ Total pressure [atm]
    %
    % --- Outputs:
    %
    % P_i ~ Partial pressure of each component [atm]

    Ftot = sum(F);

    y_i = F_i ./ Ftot;
    P_i = y_i .* Ptot;
end