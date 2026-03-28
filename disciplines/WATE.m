%% wate.m
% Weight Analysis of Turbine Engines  (WATE++ style)
%
% Sizes the engine geometry from continuity at the fan face, then estimates
% component weights using WATE++ empirical correlations scaled by disk area,
% core mass flow, and stage count.  Pylon weight is appended via a Torenbeek-
% style structural sizing that accounts for fan diameter and engine length.
%
% INPUTS  (via structs):
%   state.mdot   [kg/s]  Total engine mass flow (core + bypass), BOTH engines
%   dv.BPR       [-]     Bypass ratio
%   dv.PR_Fan    [-]     Fan pressure ratio
%   dv.PR_LPC    [-]     LPC pressure ratio
%   dv.PR_HPC    [-]     HPC pressure ratio
%   atm          struct  Atmospheric conditions at cruise
%   eng          struct  Engine reference data  (thrust_max, rho_mat)
%   wate         struct  Calibration constants  (unused directly; kept for extensibility)
%
% OUTPUTS (written into state struct):
%   state.W_engine  [kg]   Engine + pylon weight  (ONE engine)
%   state.D_fan     [m]    Fan tip diameter
%   state.D_nacelle [m]    Maximum nacelle outer diameter
%   state.L_engine  [m]    Overall engine length (fan face → nozzle exit)
%   state.A_fan     [m2]   Fan annulus area        (used by NOISE constraint)
%   state.N_stages  [-]    Total compressor stage count
%
% NOTE: Ground clearance is NOT computed here; it is evaluated in CLEARANCE.m
%       using the engine dimensions output above.
% =========================================================================

function state = wate(dv, ac, atm, eng, state, print_flag)

    if nargin < 6
        print_flag = false;
    end

    BPR    = dv.BPR;
    PR_Fan = dv.PR_Fan;
    PR_LPC = dv.PR_LPC;
    PR_HPC = dv.PR_HPC;
    OPR    = PR_Fan * PR_LPC * PR_HPC;

    % per-engine mass flow [kg/s]  (state.mdot is total for both engines)
    mdot = state.mdot / 2;
    CL = ac.CL_cr;  % use cruise CL for sizing (conservative for drag growth)
    rho = atm.rho_cr;

    % =========================================================
    % 1.  FAN FACE SIZING  (continuity, hub-to-tip ratio = 0.30)
    % =========================================================
    hub_tip_ratio = 0.30;

    % Cruise inlet total conditions (isentropic intake)
    M0   = dv.V / atm.a_cr;
    T02  = atm.T_cr * (1 + (atm.gamma - 1)/2 * M0^2);
    P02  = atm.P_cr * (T02 / atm.T_cr)^(atm.gamma / (atm.gamma - 1));

    % Static conditions at fan face (axial Mach ~ 0.55 is typical)
    M_ax   = 0.55;
    T_ax   = T02 / (1 + (atm.gamma - 1)/2 * M_ax^2);
    P_ax   = P02 / (1 + (atm.gamma - 1)/2 * M_ax^2)^(atm.gamma / (atm.gamma - 1));
    rho_ax = P_ax / (atm.R * T_ax);
    V_ax   = M_ax * sqrt(atm.gamma * atm.R * T_ax);

    % Annulus area and fan tip diameter from continuity
    A_fan  = mdot / (rho_ax * V_ax);
    D_fan  = sqrt(4 * A_fan / (pi * (1 - hub_tip_ratio^2)));

    % =========================================================
    % 2.  STAGE COUNT  (each stage ~ PR 1.25)
    % =========================================================
    psi_fan = 1.5;  % Fans are usually 1 stage, high loading
    psi_lpc = 1.25; % Low Pressure Compressor
    psi_hpc = 1.35; % High Pressure Compressor (often higher loading)

    % Calculate continuous number of stages
    n_fan = log(PR_fan) / log(psi_fan);
    n_lpc = log(PR_lpc) / log(psi_lpc);
    n_hpc = log(PR_hpc) / log(psi_hpc);
    N_stages = ceil(n_fan + n_lpc + n_hpc);  % round up to nearest whole stage
    N_comp   = N_stages - 1;                 % LPC + HPC stages

    % =========================================================
    % 3.  ENGINE LENGTH BREAKDOWN
    % =========================================================
    L_fan      = 0.50 * D_fan;              % fan module axial length [m]
    stage_pitch = 0.12;                     % axial length per compressor stage [m]
    L_combustor = 0.45;                     % combustor length [m]
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
    W_misc       = 200;   % [kg]

    W_engine_dry = W_fan + W_compressor + W_turbine + W_nacelle + W_misc;

    % =========================================================
    % 5.  PYLON WEIGHT  (Torenbeek / WATE++ structural sizing)
    % =========================================================
    n_z   = 4.5;      % ultimate load factor
    K_pyl = 0.065;    % A320-class pylon coefficient

    W_pylon_base = K_pyl * (W_engine_dry * n_z)^0.78 * (eng.thrust_max^0.12);

    % Geometry penalty: normalised to CFM56-5B reference (D=1.735 m, L=3.5 m)
    D_ref       = 1.735;
    L_ref       = 3.50;
    geom_factor = (D_fan / D_ref)^0.20 * (L_eng / L_ref)^0.30;
    W_pylon     = W_pylon_base * geom_factor;

    % =========================================================
    % 6.  TOTAL WEIGHT  (engine + pylon, per engine)
    % =========================================================
    W_engine = W_engine_dry + W_pylon;

    delta_W_eng  = ac.N_engines * (W_engine - eng.W_engine_ref);

    MTOW_new     = ac.OEW + ac.payload + ac.fuel_mass + delta_W_eng;

    S_new = MTOW_new / (0.5 * ac.CL_cr * atm.rho_cr * dv.V^2);

    W_wing = ac.W_wing* (S_new / ac.S_ref)^0.9; 

    MTOW_new = MTOW_new + W_wing - ac.W_wing;
    % =========================================================
    % 7.  WRITE OUTPUTS
    % =========================================================
    
    state.MTOW      = MTOW_new;
    state.W_engine  = W_engine;
    state.D_fan     = D_fan;
    state.D_nacelle = D_nacelle;
    state.L_engine  = L_eng;
    state.A_fan     = A_fan;
    state.N_stages  = N_stages;
    state.W_wing    = W_wing;
    state.S         = S_new;

    if print_flag
        fprintf('\n--- WATE ---\n');
        fprintf('  mdot (per eng)  = %7.2f kg/s\n', mdot);
        fprintf('  mdot_core       = %7.2f kg/s\n', mdot_core);
        fprintf('  OPR             = %7.2f\n',       OPR);
        fprintf('  N_stages        = %7d\n',         N_stages);
        fprintf('  D_fan           = %7.4f m\n',     D_fan);
        fprintf('  D_nacelle       = %7.4f m\n',     D_nacelle);
        fprintf('  L_engine        = %7.4f m\n',     L_eng);
        fprintf('  A_fan           = %7.4f m2\n',    A_fan);
        fprintf('  W_fan           = %7.1f kg\n',    W_fan);
        fprintf('  W_compressor    = %7.1f kg\n',    W_compressor);
        fprintf('  W_turbine       = %7.1f kg\n',    W_turbine);
        fprintf('  W_nacelle_skin  = %7.1f kg\n',    W_nacelle);
        fprintf('  W_engine_dry    = %7.1f kg\n',    W_engine_dry);
        fprintf('  W_pylon         = %7.1f kg\n',    W_pylon);
        fprintf('  W_engine_total  = %7.1f kg\n',    W_engine);
        fprintf('  MTOW_new        = %7.1f kg\n',    MTOW_new);
        fprintf('  W_wing          = %7.1f kg\n',    W_wing);
    end
    
end