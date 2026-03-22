%% BREGUET.m
% Breguet Range Equation
%
% Computes the aircraft cruise range using the Breguet range equation:
%
%       R = (V / (g * TSFC)) * (L/D) * ln(W_start / W_end)
%
% where:
%   W_start = MTOW  (start of cruise, after climb fuel allowance)
%   W_end   = MTOW - fuel_cruise  (end of cruise)
%
% INPUTS:
%   state.CL_CD   [-]       Lift-to-drag ratio
%   state.TSFC    [kg/N/s]  Thrust specific fuel consumption
%   state.MTOW    [kg]      Maximum take-off weight
%   dv.V          [m/s]     Cruise true airspeed
%   ac, atm, mis            Reference data
%
% OUTPUTS (written into state struct):
%   state.range   [m]       Computed cruise range
%   state.W_fuel  [kg]      Cruise fuel used
%
% =========================================================================

function state = breguet(state, dv, ac, atm, mis, print_flag)

    if nargin < 6
        print_flag = false;
    end
    V      = dv.V;
    LD     = state.CL_CD;
    TSFC   = state.TSFC;
    MTOW   = state.MTOW;
    g      = atm.g;

    % ---- 1.  FUEL FRACTIONS  ---------------------------------------------
    % Allow for take-off, climb and descent fuel burn fractions
    % (simple mission fractions, Raymer-style)
    ff_to     = 0.970;        % take-off fuel fraction
    ff_climb  = 0.985;        % climb to cruise altitude
    ff_descent= 0.990;        % descent + approach
    ff_reserve= 1 - mis.reserve_f;   % reserve fuel kept on-board

    % Weight at start of cruise
    W_start = MTOW * ff_to * ff_climb;

    % ---- 2.  BREGUET RANGE  (iterative: solve for W_end → range) ---------
    % Available cruise fuel
    W_fuel_available = ac.fuel_mass * ff_reserve - ...
                       MTOW * (1 - ff_to * ff_climb);
    W_fuel_available = max(W_fuel_available, 1);   % non-negative

    W_end   = W_start - W_fuel_available;
    W_end   = max(W_end, 0.5 * W_start);           % physical sanity bound

    % Breguet range [m]
    range   = (V / (g * TSFC)) * LD * log(W_start / W_end);

    % ---- 3.  ACCOUNT FOR DESCENT FUEL  -----------------------------------
    % Approximate final descent phase (not in Breguet) — penalise ~1%
    range   = range * ff_descent;

    % ---- 4.  WRITE OUTPUTS  ----------------------------------------------
    state.range  = range;
    state.W_fuel = W_start - W_end;

    if print_flag
        fprintf('\n--- BREGUET ---\n');
        fprintf('  RANGE       = %.1f km  (%.1f nm)\n', range/1e3, range/1852);
    end

end