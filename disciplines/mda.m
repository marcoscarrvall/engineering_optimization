function state = mda(x, data, options)
tol = options.tol;
max_iter = options.max_iter;
verbose = options.verbose;

converged = false;

% Reset state variables
state.eta_fan = 0.91;
state.D_cruise = 44103.5; 
state.objective = 0;

% Initialization 
state = thermo(state, x, data.atm, data.thermo_data, data.ac, false);
state = engine(state, x, data.atm, data.eng, data.ac, false);
state = aero(state, x, data.atm, data.eng, data.ac, false);

state.MTOW  = data.ac.OEW_ref + data.ac.W_pay + data.ac.W_fuel;  


fprintf('%-6s  %-12s  %-12s  %-12s  %-10s \n', ...
        'Iter', 'MTOW [kg]', 'rel_dMTOW', "V_cruise [m/s]", "BPR");

for iter = 1:max_iter
    MTOW_prev = state.MTOW;

    state = thermo(state, x, data.atm, data.thermo_data, data.ac, false);
    state = engine(state, x, data.atm, data.eng, data.ac, false);
    state = aero(state, x, data.atm, data.eng, data.ac, false);

    rel_change = abs(state.MTOW - MTOW_prev) / MTOW_prev;

    if verbose
        fprintf('%s\n', repmat('-', 1, 70));
        fprintf('%-6d  %-12.1f  %-10.10f  %-12.2f  %-10.2f', ...
            iter, state.MTOW, rel_change, x(1), x(2));
    end


    if rel_change < tol
        converged = true;
        break
    end

    if rel_change > 0.1
        fprintf("Warning: Divergence detected, range unchanged \n")
        state.objective = state.objective + 0.1*state.objective; 
        break
    end

end

if converged
    state = breguet(state, x, data.ac, false);
end

clearance = data.constraints.h_engine - state.D_fan/2;
state.clearance_constraint = (data.constraints.min_clearance - clearance) / data.constraints.min_clearance;
state.TIT_constraint = (data.constraints.TIT_max - state.TIT) / data.constraints.TIT_max;
state.M_tip_constraint = (data.constraints.m_tip_max - state.M_tip) / data.constraints.m_tip_max;