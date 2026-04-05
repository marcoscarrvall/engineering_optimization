function state = breguet(state, x, atm, ac, mission, print_flag)

    if nargin < 6
        print_flag = false;
    end
    V      = x(1);
    LD     = state.CL_CD;
    TSFC   = state.TSFC;
    MTOW   = state.MTOW;
    g      = atm.g;
    fprintf('\n  Breguet: V=%.1f m/s, L/D=%.2f, TSFC=%.4f, MTOW=%.0f => ', V, LD, TSFC*10000, MTOW);

    % ---- 1.  FUEL FRACTIONS  ---------------------------------------------
    % Allow for take-off, climb and descent fuel burn fractions
    % (simple mission fractions, Raymer-style)
    ff_to = mission.ff_takeoff;           % take-off fuel fraction
    ff_climb = mission.ff_climb;     % climb fuel fraction
    ff_descent = mission.ff_descent; % descent fuel fraction

    % Weight at start of cruise
    W_start = MTOW * ff_to * ff_climb;

    W_end   = W_start - ac.W_fuel;

    range   = (V / (g * TSFC)) * LD * log(W_start / W_end);

    range   = range * ff_descent;

    state.range  = range;

    if print_flag
        fprintf('\n--- BREGUET ---\n');
        fprintf('  RANGE       = %.1f km  (%.1f nm)\n', range/1e3, range/1852);
    end

end