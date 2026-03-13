function [Emissions] = Emissions(Pt3, Tt3, f, tau_res, Phi, m_dot_fuel)
    % Connection: 4: P_{T3}, T_{T3}, f, \tau_{res}, \Phi
    
    % Constants and Reference Values
    P_ref = 101325; % [Pa]
    T_ref = 288.15; % [K]
    
    %% 1. NOx Calculation (Schwartz Dallara Method)
    % a, b, m are empirical constants for a modern lean-burn combustor
    a = 0.035; 
    b = 288;
    m = 1.4;
    
    EI_NOx = a * (Pt3/P_ref)^0.5 * exp((Tt3 - T_ref)/b) * (Phi^m)* tau_res;
    
    %% 2. Soot Calculation
    % Soot index increases with pressure (Pt3) 
    % and is scaled by the fuel-air ratio (f)
    K_soot = 0.001; % Empirical scaling
    EI_Soot = K_soot * (Pt3/P_ref)^1.2 * f; 
    
    %% 3. CO2 Calculation
    % Strictly dependent on fuel chemistry
    EI_CO2 = 3150; % g/kg of fuel
    
    %% 4. Total Mass Flow of Pollutants [g/s]
    % m_dot_fuel must be in kg/s
    Emissions.NOx_rate  = EI_NOx * m_dot_fuel;
    Emissions.Soot_rate = EI_Soot * m_dot_fuel;
    Emissions.CO2_rate  = EI_CO2 * m_dot_fuel;

end