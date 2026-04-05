function [c_noise] = noise(state, x, atm, thermo, eng)

    Cp_c  = thermo.Cp_air;
    gam_c = thermo.gamma_c;
    eta_f = state.eta_fan;
    PR_fan = x(3);

    % ---- Fan inlet total conditions ------------------------------------
    % Ram temperature rise from flight velocity
    a_cr = sqrt(atm.gamma * atm.R * atm.T_cruise);
    M0   = x(1) / a_cr;                                        % flight Mach
    T02  = atm.T_cruise * (1 + (gam_c - 1)/2 * M0^2);         % total temp at fan inlet

    % ---- Fan isentropic temperature rise from PR_fan (x3) --------------
    T021_is = T02 * x(3)^((gam_c - 1) / gam_c);               % ideal outlet total temp
    T021    = T02 + (T021_is - T02) / eta_f;                   % actual outlet total temp

    % ---- Fan specific work scaled by BPR (x2) --------------------------
    % For fixed thrust, a higher BPR engine moves more air at lower
    % velocity — the fan does less work per unit of TOTAL mass flow.
    % This physically links BPR to tip speed and hence noise.
    BPR     = x(2);
    delta_h = Cp_c * (T021 - T02) / (1 + BPR);                % [J/kg] per total flow

    % ---- Tip speed from Euler work equation ----------------------------
    psi_fan = 0.40;                                             % work coefficient [-]
    U_tip   = sqrt(delta_h / psi_fan);                         % [m/s]

    % ---- Static conditions at fan face (axial Mach = 0.55) -------------
    M_ax  = 0.55;
    T_ax  = T02  / (1 + (gam_c - 1)/2 * M_ax^2);              % static temp at fan face
    a_ax  = sqrt(gam_c * atm.R * T_ax);                        % static speed of sound
    V_ax  = M_ax * a_ax;                                       % axial velocity [m/s]

    % ---- Relative tip Mach number --------------------------------------
    % V_rel combines rotational tip speed and axial velocity.
    % Referenced to STATIC speed of sound (physically correct for blade noise).
    V_rel = sqrt(U_tip^2 + V_ax^2);
    M_tip = V_rel / a_ax;
    M_tip = state.M_tip;  % override with stored value for consistent reporting

    % ---- Store and evaluate constraint ---------------------------------
    state.M_tip = M_tip;

    % Normalised inequality: g(x) = (M_tip - M_tip_max) / M_tip_max <= 0
    c_noise  = (M_tip - eng.M_tip_max) / eng.M_tip_max;
    violated = c_noise > 0;
    margin   = -c_noise;       % positive = safe headroom

    if violated
        fprintf('  STATUS        : *** VIOLATED (fan too fast / too loud) ***\n');
    end

end