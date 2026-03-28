clear; clc;

data = TestAC_data;
x.BPR = 14.0;
x.V   = sqrt(data.atm.gamma * data.atm.R * data.atm.T_cruise) * 0.78;

x_consts.PR_fan = 1.6;
x_consts.PR_LPC = 2.6;
x_consts.PR_HPC = 6.1;

atm         = data.atm;
thermo_data = data.thermo;
eng         = data.eng;
wate_consts = data.wate;
ac          = data.ac;

%% ---- Initialisation ----------------------------------------------------
state        = data.state;
ac.S_ref  = (ac.MTOW - ac.W_fuel/2) * atm.g / ...
            (0.5 * ac.CL_cruise * atm.rho_cruise * x.V^2);
state.D_fan  = eng.D_fan_ref;
state.L_eng  = eng.L_eng_ref;

% Approximate initial thrust from MTOW and assumed L/D = 17.5
state.D_cruise  = ac.MTOW * atm.g / 17.5;

% Approximate initial mdot from thermo at that thrust
state.MTOW      = ac.MTOW;
state           = thermo(state, x, x_consts, atm, thermo_data, ac);

%% ---- MDA loop ----------------------------------------------------------
tol        = 1e-4;   
max_iter   = 50;
converged  = false;

fprintf('%-6s  %-12s  %-12s  %-12s  %-10s  %-10s\n', ...
    'Iter', 'MTOW [kg]', 'D_cruise [N]', 'mdot [kg/s]', 'D_fan [m]', 'rel_dMTOW');
fprintf('%s\n', repmat('-', 1, 70));

for iter = 1:max_iter
    MTOW_prev = state.MTOW;

    state = wate  (state, x, x_consts, atm, eng, wate_consts, ac);
    state = aero  (state, x, atm, eng, ac);
    state = thermo(state, x, x_consts, atm, thermo_data, ac);

    rel_change = abs(state.MTOW - MTOW_prev) / MTOW_prev;

    fprintf('%-6d  %-12.1f  %-12.1f  %-12.2f  %-10.4f  %-10.2e\n', ...
        iter, state.MTOW, state.D_cruise, state.mdot, state.D_fan, rel_change);

    if rel_change < tol
        converged = true;
        break
    end
end

%% ---- Result ------------------------------------------------------------
fprintf('%s\n', repmat('-', 1, 70));
if converged
    fprintf('Converged in %d iterations (rel_dMTOW < %.0e)\n\n', iter, tol);
else
    fprintf('WARNING: did not converge in %d iterations\n\n', max_iter);
end
