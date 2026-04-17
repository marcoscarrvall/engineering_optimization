function state = engine(state, x, atm, eng, ac, print_flag)

    if nargin < 6
        print_flag = false;
    end

    BPR    = x(2);
    PR_fan = x(3);
    PR_LPC = x(4);
    PR_HPC = x(5);

    mdot = state.mdot_total / 2;
    atm.gamma = 1.4;

    M0   = x(1) / sqrt(atm.gamma * atm.R * atm.T_cruise);
    T02  = atm.T_cruise * (1 + (atm.gamma - 1)/2 * M0^2);
    P02  = atm.P_cruise * (T02 / atm.T_cruise)^(atm.gamma / (atm.gamma - 1));

    M_ax   = 0.5; % Assume flow is reduced to M=0.5 at fan face
    T_ax   = T02 / (1 + (atm.gamma - 1)/2 * M_ax^2);
    P_ax   = P02 / (1 + (atm.gamma - 1)/2 * M_ax^2)^(atm.gamma / (atm.gamma - 1));
    rho_ax = P_ax / (atm.R * T_ax);
    V_ax   = M_ax * sqrt(atm.gamma * atm.R * T_ax);

    % Annulus area and fan tip diameter from continuity
    K = 2.15*BPR^0.3;  % Correction factor due to flow non-uniformity and blockage
    A_fan  = K * mdot / (rho_ax * V_ax);
    D_fan  = sqrt(4 * A_fan / pi);

    % Calculate continuous number of stages
    n_fan = 1;
    n_lpc = PR_LPC / log(2.5);
    n_hpc = PR_HPC / log(46.4);
    N_stages = n_fan + n_lpc + n_hpc;  % We don't round this - physically wrong but acceptable for this model
    N_comp   = N_stages - 1;                 % LPC + HPC stages

    % 3.  ENGINE LENGTH BREAKDOWN
    L_fan      = 0.75 * D_fan;              % fan module axial length [m]
    stage_pitch = 0.3;                     % axial length per compressor stage [m]
    L_combustor = 0.75;                     % combustor length [m]
    L_core      = N_comp * stage_pitch + L_combustor;
    L_nozzle    = 0.6 * L_core;            % nozzle / mixer section [m]
    L_eng       = 3/5 * (L_fan + L_core + L_nozzle); % Correction factor 3/5

    % Weight Components
    W_fan        = 0.07 * (D_fan^0.1) * eng.rho_mat * (1 + 0.10 * BPR^1.0);
    W_compressor = 103.61 * (N_comp * 0.65);
    W_turbine    = 153.92 * (N_comp * 0.35);
    W_nacelle    = 164.27 * pi * D_fan^0.1;
    W_misc       = 950;  

    W_engine = W_fan + W_compressor + W_turbine + W_nacelle + W_misc;

    delta_W_eng  = ac.N_eng * (W_engine - eng.W_eng_ref);

    MTOW_new     = (ac.OEW_ref + ac.W_pay + ac.W_fuel + delta_W_eng);

    M_tip = 0.55 * sqrt((1+BPR)*(PR_fan-1));
    
    state.MTOW      = MTOW_new;
    state.W_engine  = W_engine;
    state.D_fan     = D_fan;
    state.L_eng     = L_eng;
    state.A_fan     = A_fan;
    state.N_stages  = N_stages;
    state.M_tip     = M_tip;
    state.eta_fan   = eng.eta_fan-0.1*(max(0, M_tip-1.2))^2; % simple model: fan efficiency drops if we have shock waves

    if print_flag
        fprintf('\n--- WATE ---');
        fprintf("\n Fan diameter (D_fan) = %5.3f m", D_fan);
        fprintf("\n Engine length (L_eng) = %5.3f m", L_eng);
        fprintf("\n Engine weight (W_engine) = %5.1f kg", W_engine);
        fprintf("\n Change in engine weight (delta_W_eng) = %5.1f kg", delta_W_eng);
        fprintf("\n Fan tip Mach number (M_tip) = %5.3f", M_tip);
        fprintf("\n Fan isentropic efficiency (eta_fan) = %5.3f\n", state.eta_fan);
    end
    
end