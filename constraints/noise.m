%% NOISE.m
% Fan Tip Mach Number  —  Noise Constraint
%
% Computes the fan-tip RELATIVE Mach number (M_tip) from the engine geometry
% (WATE) and fan thermodynamic operating point, then checks it against the
% allowable limit.
%
% Physics:
%   The fan-tip speed U_tip follows from the Euler turbomachinery equation:
%       delta_h_fan = U_tip^2 * psi_fan
%   where psi_fan is the fan work coefficient (typically 0.35–0.45).
%
%   The relative tip Mach number combines blade speed and axial velocity:
%       V_rel = sqrt(U_tip^2 + V_ax^2)
%       M_tip = V_rel / a_inlet
%
% Constraint is NORMALISED by M_tip_max so the value is dimensionless
% and O(1) for use in gradient-based optimisers (e.g. fmincon):
%
%   c_noise = (M_tip - M_tip_max) / M_tip_max  <=  0
%
% A value of +0.05 means M_tip is 5% above the limit;
% a value of -0.10 means there is 10% headroom below the limit.
%
% INPUTS:
%   state.D_fan   [m]    Fan tip diameter (from WATE)
%   state.A_fan   [m2]   Fan annulus area (from WATE)
%   dv.PR_Fan     [-]    Fan pressure ratio
%   dv.V          [m/s]  Cruise velocity
%   atm, thermo   structs
%   con.M_tip_max [-]    Maximum allowable tip Mach number
%
% OUTPUTS:
%   c_noise   [-]   Normalised constraint  (<=0 satisfied, >0 violated)
%   violated  [bool]
%   margin    [-]   Normalised headroom  = -c_noise  (positive = safe)
%   state.M_tip [-]  Written back for diagnostics and post-processing
%
% =========================================================================

function [c_noise] = noise(state, dv, atm, thermo, con)

    Cp_c  = thermo.Cp_air;
    gam_c = thermo.gamma_c;
    eta_f = thermo.eta_fan;

    % ---- Fan inlet total conditions (same as THERMO station 2) ----------
    M0   = dv.V / atm.a_cr;
    T02  = atm.T_cr * (1 + (gam_c - 1)/2 * M0^2);

    % ---- Fan temperature rise from pressure ratio (polytropic) ----------
    T021_is = T02 * dv.PR_Fan^((gam_c - 1) / gam_c);
    T021    = T02 + (T021_is - T02) / eta_f;
    delta_h = Cp_c * (T021 - T02);              % [J/kg] fan specific work

    % ---- Tip speed from Euler work equation -----------------------------
    psi_fan = 0.40;
    U_tip   = sqrt(delta_h / psi_fan);          % [m/s]

    % ---- Axial velocity at fan face (from continuity, same as WATE) -----
    M_ax = 0.55;
    T_ax = T02 / (1 + (gam_c - 1)/2 * M_ax^2);
    V_ax = M_ax * sqrt(gam_c * atm.R * T_ax);

    % ---- Relative tip velocity and Mach number --------------------------
    V_rel  = sqrt(U_tip^2 + V_ax^2);
    a_in   = sqrt(gam_c * atm.R * T02);         % speed of sound at fan inlet total
    M_tip  = V_rel / a_in;

    state.M_tip = M_tip;

    % Normalised inequality: g(x) = (M_tip - M_tip_max) / M_tip_max <= 0
    c_noise  = (M_tip - con.M_tip_max) / con.M_tip_max;
    violated = c_noise > 0;
    margin   = -c_noise;       % [-]  positive = safe headroom


    if violated
        fprintf('  STATUS        : *** VIOLATED (fan too fast / too loud) ***\n');
    else
        fprintf('  STATUS        : Satisfied\n');
        fprintf('  Margin        : %+.6f [-]\n', margin);
    end

end