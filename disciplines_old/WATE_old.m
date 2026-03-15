function [Total_Weight,delta_CD,delta_CL] = WATE_old(m_dot_total, BPR, OPR, rho_mat, thrust_max, rho, V, S_ref, AR, e)

    % Connection: 3: \dot{m}_{total}, BPR, D_{fan}, N_{stages}, \rho_{mat}
    hub_tip_ratio = 0.3; % Typical for high-bypass fans
    Area_required = m_dot_total / (rho * V);
    D_fan = sqrt( (4 * Area_required) / (pi * (1 - hub_tip_ratio^2)) );

    %% 1. Component Mass Flow Split
    m_dot_core = m_dot_total / (1 + BPR);

    N_stages = ceil(log(OPR) / log(1.25));
    %% 2. WATE++ Component Breakdown
    
    % A. Fan Module: Scaled by Disk Area and BPR
    % Based on WATE++ correlation for large high-bypass fans
    W_fan = 0.085 * (D_fan^2.1) * rho_mat * (1 + 0.1 * BPR);
    
    % B. Core Compressor (LPC + HPC): Scaled by stage count and core flow
    % WATE++ assumes weight increases with the number of stages and air density
    W_compressor = 1.15 * m_dot_core * (N_stages * 0.65);
    
    % C. Turbine & Combustor Module:
    % Heavily influenced by the high-pressure core mass flow
    W_turbine = 0.95 * m_dot_core * (N_stages * 0.35);
    
    % D. Nacelle & Accessories:
    % WATE++ uses a scaling factor based on the fan diameter (D_fan)
    W_nacelle = 0.045 * pi * D_fan^2 * 150; % Structural weight per m^2
    W_misc = 200; % Control systems, pumps, and fluids
    
    
    %% 3. Summation
    Total_Weight_engine = W_fan + W_compressor + W_turbine + W_nacelle + W_misc;
    
    % n_z: Load factor
    n_z = 4.5; 
    
    % K factor for A320 class
    K_pyl = 0.065; 

    % 1. Baseline Weight (Your existing formula)
    W_base = K_pyl * (Total_Weight_engine * n_z)^0.78 * (thrust_max^0.12);

    % 2. Geometry Penalty (Characteristic Adjustment)
    % A wider fan (D_fan) increases frontal area/drag loads.
    % A longer engine (L_engine) increases the bending moment.
    % We use a non-dimensional ratio compared to a reference A320 engine (e.g., D=1.7m, L=3.5m)
    L_fan = 0.5 * D_fan;

    % 2. Core Length
    % Each stage typically takes up 0.1 to 0.15 meters of axial space 
    % plus the combustor length (approx 0.4m)
    stage_pitch = 0.12; % [m]
    L_combustor = 0.45; % [m]
    L_core = (N_stages * stage_pitch) + L_combustor;

    % 3. Nozzle & Integration Factor
    % High BPR engines have shorter, fatter nacelles
    L_nozzle = 0.25 * L_core;

    % 4. Total Length
    L_eng = L_fan + L_core + L_nozzle;

    geom_factor = (D_fan / 1.7)^0.2 * (L_eng / 3.5)^0.3;

    % 3. Final Pylon Weight
    W_pylon = W_base * geom_factor;

    % 4. Total Engine + Pylon Weight
    Total_Weight = Total_Weight_engine + W_pylon;

    S_wet_ratio = (D_fan * L_eng) / (D_ref * L_ref);

    delta_CD0 = 0.0020 * (S_wet_ratio - 1); % Example increment
    
    % 2. Interference Drag (increases rapidly as diameter grows)
    delta_CD_int = 0.0005 * (D_fan / D_ref)^2.5; 
    
    delta_CD_parasitic = delta_CD0 + delta_CD_int;

    delta_CL = (Total_Weight - Total_Weight_ref) / (0.5 * rho * V^2 * S_ref);

    delta_CD_induced = delta_CL^2 / (pi * AR * e);

    delta_CD = delta_CD_parasitic + delta_CD_induced;


end