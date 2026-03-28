function state = aero(state, x, atm, eng, ac, print_flag)

    if nargin < 6
        print_flag = false;
    end

    V    = x.V;
    rho  = atm.rho_cruise;
    q    = 0.5 * rho * V^2;          % dynamic pressure [Pa]

    D_ref = eng.D_fan_ref;
    L_ref = eng.L_eng_ref;
    CD_total = ac.CD_cruise;         % baseline total drag coefficient (from data)

    % Wetted-area ratio of current nacelle vs. reference
    S_wet_ratio   = 2*(state.D_fan * state.L_eng-D_ref * L_ref) / (ac.S_ref);

    % Zero-lift drag increment from changed nacelle size
    delta_CD_parasitic = ac.CD_parasitic_ref *S_wet_ratio;

    % Interference drag — scales strongly with fan diameter growth
    delta_CD_int  = 0.0005 * (state.D_fan / D_ref)^2.5;

    delta_CD0 = delta_CD_parasitic + delta_CD_int;
    D_total  = (CD_total + delta_CD0) * q * state.S;                % total from scratch

    state.D_cruise = D_total;
    state.CL_CD   = ac.CL_cruise / (CD_total+delta_CD0);

    if print_flag
        fprintf('\n--- AERO ---\n');
        fprintf('  L/D             = %8.4f\n', state.CL_CD);
        fprintf('  Drag (total)    = %8.1f N\n', D_total);
    end

end