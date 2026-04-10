function state = aero(state, x, atm, eng, ac, print_flag)

    if nargin < 6
        print_flag = false;
    end

    V    = x(1);
    rho  = atm.rho_cruise;
    q    = 0.5 * rho * V^2;          % dynamic pressure [Pa]

    D_ref = eng.D_fan_ref;
    L_ref = eng.L_eng_ref;
    CD_total_ref = ac.CD_cruise_ref;         % baseline total drag coefficient (from data)

    S_wet_ratio   = (state.D_fan * state.L_eng-D_ref * L_ref + ac.S_ref) / (ac.S_ref);
    delta_CD_parasitic = ac.CD_parasitic_ref * (S_wet_ratio-1);
    delta_CD_int  = 2 * (state.D_fan / D_ref - 1)^2.5;
    delta_CD0 = (delta_CD_parasitic + delta_CD_int) * (V/ 230); % Scale profile drag with velocity 
    if V > 250 
        delta_CDwave = 0.0005 * (V-250)/10; 
    else 
        delta_CDwave = 0; 
    end
    delta_CD0 = delta_CD0 + delta_CDwave;

    CL = state.MTOW*9.81 / (q * ac.S_ref);
    delta_CL2 = CL^2 - ac.CL_cruise_ref^2;
    delta_CDi = delta_CL2/ (pi * ac.AR * ac.e);

    delta_CD = delta_CD0 + delta_CDi;
    CD = CD_total_ref + delta_CD;

    D_total  = CD * q * ac.S_ref;                

    state.D_cruise = D_total;
    state.CL_CD   = CL / CD;

    if print_flag
        fprintf('\n--- AERO ---');
        fprintf("\n Lift-to-drag ratio (CL/CD) = %5.1f\n", state.CL_CD);
    end

end