%% BUILD_REFERENCE.m
% Shared helper for all test files.
%
% Returns the same reference structs that data.m populates, but as a
% function so each test file can call it without relying on workspace state.
%
% Usage:
%   [atm, dv, ac, con, eng, thermo_data, wate_data, mis, state] = build_reference();
% =========================================================================

function [atm, dv, ac, con, eng, thermo_data, wate_data, mis, state] = build_reference()

    % ---- Atmosphere (ISA, 35 000 ft) ------------------------------------
    atm.T_sl   = 288.15;
    atm.P_sl   = 101325;
    atm.rho_sl = 1.225;
    atm.gamma  = 1.4;
    atm.R      = 287.05;
    atm.g      = 9.80665;
    atm.h_cr   = 10668;
    atm.T_cr   = 218.81;
    atm.P_cr   = 23842;
    atm.rho_cr = 0.3796;
    atm.a_cr   = sqrt(atm.gamma * atm.R * atm.T_cr);

    % ---- Design vector --------------------------------------------------
    dv.V      = 230;
    dv.BPR    = 5.9;
    dv.PR_HPC = 9.5;
    dv.PR_LPC = 2.6;
    dv.PR_Fan = 1.7;
    dv.OPR    = dv.PR_Fan * dv.PR_LPC * dv.PR_HPC;

    % ---- Aircraft (A320) ------------------------------------------------
    ac.MTOW       = 77000;
    ac.OEW        = 42600;
    ac.fuel_mass  = 20100;
    ac.payload    = 14300;
    ac.S_ref      = 122.4;
    ac.b          = 34.1;
    ac.AR         = ac.b^2 / ac.S_ref;
    ac.e          = 0.82;
    ac.sweep_LE   = 27.1;
    ac.CL_cr      = 0.56;
    ac.CD0_clean  = 0.0220;
    ac.CL_CD_ref  = 17.5;
    ac.N_engines  = 2;

    % ---- Engine reference (CFM56-5B) ------------------------------------
    eng.N_engines    = 2;
    eng.mdot_ref     = 350;
    eng.T_sl_ref     = 120000;
    eng.thrust_max   = 120000;
    eng.TSFC_ref     = 1.65e-5;
    eng.W_engine_ref = 2350;
    eng.D_fan_ref    = 1.735;
    eng.L_engine_ref = 3.50;
    eng.rho_mat      = 4430;
    eng.pylon_h      = 0.40;
    eng.wing_h_AGL   = 1.80;

    % ---- Thermodynamic constants ----------------------------------------
    thermo_data.Cp_air     = 1005;
    thermo_data.Cp_gas     = 1150;
    thermo_data.gamma_c    = 1.40;
    thermo_data.gamma_h    = 1.33;
    thermo_data.LHV        = 43.2e6;
    thermo_data.eta_cc     = 0.995;
    thermo_data.dP_cc_frac = 0.05;
    thermo_data.eta_fan    = 0.92;
    thermo_data.eta_LPC    = 0.90;
    thermo_data.eta_HPC    = 0.87;
    thermo_data.eta_HPT    = 0.91;
    thermo_data.eta_LPT    = 0.93;
    thermo_data.eta_mech   = 0.99;
    thermo_data.FAR        = 0.0270;

    % ---- Constraints ----------------------------------------------------
    con.TIT_max       = 1700;
    con.clearance_min = 0.45;
    con.M_tip_max     = 1.50;

    % ---- WATE calibration -----------------------------------------------
    wate_data.placeholder = 0;

    % ---- Mission profile ------------------------------------------------
    mis.R_design  = 3300e3;
    mis.h_cr      = atm.h_cr;
    mis.reserve_f = 0.05;

    % ---- Initial MDA state ----------------------------------------------
    state.mdot      = eng.mdot_ref;
    state.W_engine  = eng.W_engine_ref;
    state.MTOW      = ac.MTOW;
    state.D_total   = 60000 * 2;
    state.CL_CD     = ac.CL_CD_ref;
    state.TSFC      = eng.TSFC_ref;
    state.D_fan     = eng.D_fan_ref;
    state.D_nacelle = eng.D_fan_ref * 1.10;
    state.L_engine  = eng.L_engine_ref;
    state.A_fan     = 0;
    state.N_stages  = 0;
    state.TIT       = 0;
    state.range     = 0;
    state.M_tip     = 0;
    state.CL        = 0;
    state.clearance = 0;

end