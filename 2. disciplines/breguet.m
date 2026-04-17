function state = breguet(state, x, ac, print_flag)

    if nargin < 4
        print_flag = false;
    end

    V      = x(1);

    range   = (V / (9.81 * state.TSFC)) * state.CL_CD * log(state.MTOW / (state.MTOW - ac.W_fuel));

    state.range  = range;

    state.objective = -(range-ac.range)/ac.range;

    if print_flag
        fprintf('\n--- BREGUET ---');
        fprintf('\n  RANGE       = %.1f km', range/1e3);
        fprintf('\n  Objective    = %5.3f\n', state.objective);
    end

end