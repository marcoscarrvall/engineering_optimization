function state = wate(state, x, x_consts, atm, eng, wate_consts, ac, print_flag)

    if nargin < 8
        print_flag = false;
    end

    BPR    = x(2);
    PR_fan = x_consts.PR_fan;
    PR_LPC = x_consts.PR_LPC;
    PR_HPC = x_consts.PR_HPC;
    OPR    = PR_fan * PR_LPC * PR_HPC;

    % per-engine mass flow [kg/s]  (state.mdot is total for both engines)
    mdot = state.mdot / 2*0.7;

    hub_tip_ratio = wate_consts.hub_to_tip;  % typical turbofan value

    % Cruise inlet total conditions (isentropic intake)
    M0   = x(1) / sqrt(atm.gamma * atm.R * atm.T_cruise);
    T02  = atm.T_cruise * (1 + (atm.gamma - 1)/2 * M0^2);
    P02  = atm.P_cruise * (T02 / atm.T_cruise)^(atm.gamma / (atm.gamma - 1));

    % Static conditions at fan face (axial Mach ~ 0.55 is typical)
    M_ax   = wate_consts.M_axial;
    T_ax   = T02 / (1 + (atm.gamma - 1)/2 * M_ax^2);
    P_ax   = P02 / (1 + (atm.gamma - 1)/2 * M_ax^2)^(atm.gamma / (atm.gamma - 1));
    rho_ax = P_ax / (atm.R * T_ax);
    V_ax   = M_ax * sqrt(atm.gamma * atm.R * T_ax);

    % Annulus area and fan tip diameter from continuity
    A_fan  = mdot / (rho_ax * V_ax);
    D_fan  = sqrt(4 * A_fan / (pi * (1 - hub_tip_ratio^2)));

    psi_lpc = 1.25; % Low Pressure Compressor
    psi_hpc = 1.35; % High Pressure Compressor (often higher loading)

    % Calculate continuous number of stages
    n_fan = 1;
    n_lpc = log(PR_LPC) / log(psi_lpc);
    n_hpc = log(PR_HPC) / log(psi_hpc);
    N_stages = ceil(n_fan + n_lpc + n_hpc);  % round up to nearest whole stage
    N_comp   = N_stages - 1;                 % LPC + HPC stages

    % =========================================================
    % 3.  ENGINE LENGTH BREAKDOWN
    % =========================================================
    L_fan      = 0.50 * D_fan;              % fan module axial length [m]
    stage_pitch = wate_consts.p_stage;                     % axial length per compressor stage [m]
    L_combustor = wate_consts.L_combustor;                     % combustor length [m]
    L_core      = N_comp * stage_pitch + L_combustor;
    L_nozzle    = 0.25 * L_core;            % nozzle / mixer section [m]
    L_eng       = L_fan + L_core + L_nozzle;

    % Nacelle outer diameter (10 % clearance over fan tip)
    D_nacelle = D_fan * 1.10;

    % =========================================================
    % 4.  COMPONENT WEIGHTS  (WATE++ correlations)
    % =========================================================
    mdot_core = mdot / (1 + BPR);

    % Fan — disk area and BPR loading
    W_fan        = 0.085 * (D_fan^2.1) * eng.rho_mat * (1 + 0.10 * BPR);

    % Compressor (LPC + HPC) — stage count and core flow
    W_compressor = 1.15 * mdot_core * (N_comp * 0.65);

    % Turbine + combustor — lighter per stage than compressor
    W_turbine    = 0.95 * mdot_core * (N_comp * 0.35);

    % Nacelle skin — wetted area × structural areal density (150 kg/m²)
    W_nacelle    = 0.045 * pi * D_fan^2 * 150;

    % Accessories — controls, pumps, fluids (fixed allowance)
    W_misc       = wate_consts.W_misc;   % [kg]

    W_engine = W_fan + W_compressor + W_turbine + W_nacelle + W_misc;

    delta_W_eng  = ac.N_eng * (W_engine - eng.W_eng_ref);

    MTOW_new     = ac.OEW + ac.W_pay + ac.W_fuel + delta_W_eng;


    S_new = (MTOW_new-ac.W_fuel/2)*9.81 / (0.5 * ac.CL_cruise * atm.rho_cruise * x(1)^2);


    W_wing = ac.W_wing* (S_new / ac.S_ref)^0.9; 

    MTOW_new = MTOW_new + W_wing - ac.W_wing;
    % =========================================================
    % 7.  WRITE OUTPUTS
    % =========================================================
    
    state.MTOW      = MTOW_new;
    state.W_engine  = W_engine;
    state.D_fan     = D_fan;
    state.D_nacelle = D_nacelle;
    state.L_eng  = L_eng;
    state.A_fan     = A_fan;
    state.N_stages  = N_stages;
    state.W_wing    = W_wing;
    state.S         = S_new;

    if print_flag
        fprintf('\n--- WATE ---\n');
        fprintf('  N_stages        = %7d\n',         N_stages);
        fprintf('  D_fan           = %7.4f m\n',     D_fan);
        fprintf('  L_eng           = %7.4f m\n',     L_eng);
        fprintf('  W_engine_total  = %7.1f kg\n',    W_engine);
        fprintf('  MTOW_new        = %7.1f kg\n',    MTOW_new);
    end
    
end