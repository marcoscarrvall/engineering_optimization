function [c_noise] = noise(state, x, x_consts, atm, thermo, eng)

    Cp_c  = thermo.Cp_air;
    gam_c = thermo.gamma_c;
    eta_f = thermo.eta_fan;

    % ---- Fan inlet total conditions (same as THERMO station 2) ----------
    a_cr = sqrt(atm.gamma * atm.R * atm.T_cruise);
    M0   = x(1) / a_cr;
    T02  = atm.T_cruise * (1 + (gam_c - 1)/2 * M0^2);

    % ---- Fan temperature rise from pressure ratio (polytropic) ----------
    T021_is = T02 * x_consts.PR_fan^((gam_c - 1) / gam_c);
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
    c_noise  = (M_tip - eng.M_tip_max) / eng.M_tip_max;
    violated = c_noise > 0;
    margin   = -c_noise;       % [-]  positive = safe headroom


    if violated
        fprintf('  STATUS        : *** VIOLATED (fan too fast / too loud) ***\n');
    else
        fprintf('  STATUS        : Satisfied\n');
        fprintf('  Margin        : %+.6f [-]\n', margin);
    end

end