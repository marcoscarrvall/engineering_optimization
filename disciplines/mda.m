function state = mda(x, x_consts, data, options)
tol = options.tol;
max_iter = options.max_iter;
verbose = options.verbose;

%% ---- 1. Unpack data ----------------------------------------------------
atm         = data.atm;
thermo_data = data.thermo;
eng         = data.eng;
wate_consts = data.wate;
ac          = data.ac;

% Override S_ref to be consistent with cruise lift balance at mid-cruise weight
ac.S_ref = (ac.MTOW - ac.W_fuel/2) * atm.g / ...
           (0.5 * ac.CL_cruise * atm.rho_cruise * x(1)^2);

%% ---- 2. Initialisation -------------------------------------------------
state           = data.state;
state.MTOW      = ac.MTOW;
state.D_fan     = eng.D_fan_ref;
state.L_eng     = eng.L_eng_ref;
state.S         = ac.S_ref;

% Initial drag from MTOW and assumed L/D
state.D_cruise  = ac.MTOW * atm.g / 17.5;

% Initial mdot from thermo at that thrust
state = thermo(state, x, x_consts, atm, thermo_data, ac);

%% ---- 3. MDA loop -------------------------------------------------------
if verbose
    fprintf('%-6s  %-12s  %-12s  %-12s  %-10s  %-10s\n', ...
        'Iter', 'MTOW [kg]', 'D_cruise [N]', 'mdot [kg/s]', 'D_fan [m]', 'rel_dMTOW');
    fprintf('%s\n', repmat('-', 1, 70));
end

converged = false;

for iter = 1:max_iter
    MTOW_prev = state.MTOW;

    state = wate  (state, x, x_consts, atm, eng, wate_consts, ac);
    state = aero  (state, x, atm, eng, ac);
    state = thermo(state, x, x_consts, atm, thermo_data, ac);

    rel_change = abs(state.MTOW - MTOW_prev) / MTOW_prev;

    if verbose
        fprintf('%-6d  %-12.1f  %-12.1f  %-12.2f  %-10.4f  %-10.2e\n', ...
            iter, state.MTOW, state.D_cruise, state.mdot, state.D_fan, rel_change);
    end

    if rel_change < tol
        converged = true;
        break
    end
end

%% ---- 4. Convergence report ---------------------------------------------
if verbose
    fprintf('%s\n', repmat('-', 1, 70));
    if converged
        fprintf('Converged in %d iterations (rel_dMTOW < %.0e)\n\n', iter, tol);
    else
        fprintf('WARNING: did not converge in %d iterations\n\n', max_iter);
    end
end

if ~converged
    warning('mda:noConvergence', ...
        'MDA did not converge in %d iterations. Results may be unreliable.', max_iter);
end

end