%% TEST_BREGUET.m
% Unit tests for BREGUET.m
%
% Tests that the Breguet range equation produces physically consistent
% cruise range and fuel-burn outputs.
% =========================================================================

function test_breguet()
    fprintf('========================================\n');
    fprintf(' TEST SUITE: BREGUET.m\n');
    fprintf('========================================\n');
    pass = 0; fail = 0;

    % ---- Build reference structs ----------------------------------------
    [atm, dv, ac, ~, ~, ~, ~, mis, state] = build_reference();

    % Provide analysis outputs (normally from THERMO and AERO)
    state.CL_CD = 17.5;           % reference L/D
    state.TSFC  = 1.65e-5;        % [kg/N/s] reference TSFC
    state.MTOW  = ac.MTOW;        % [kg]

    % ---- Run module -------------------------------------------------------
    state = breguet(state, dv, ac, atm, mis);

    % ---- Test 1: Range positive ----------------------------------------
    [pass, fail] = check('range > 0', state.range > 0, pass, fail);

    % ---- Test 2: Range in realistic bounds for A320 --------------------
    % A320 typical range ~2700–6000 km depending on payload and fuel
    [pass, fail] = check('range in [500 km, 10000 km]', ...
        state.range >= 500e3 && state.range <= 10000e3, pass, fail);

    % ---- Test 3: Range within factor of 2 of design range --------------
    [pass, fail] = check('range within factor-2 of design range (3300 km)', ...
        state.range >= mis.R_design / 2 && state.range <= mis.R_design * 2, pass, fail);

    % ---- Test 4: Fuel used positive and less than available fuel --------
    [pass, fail] = check('W_fuel > 0', state.W_fuel > 0, pass, fail);
    [pass, fail] = check('W_fuel < total fuel mass', ...
        state.W_fuel < ac.fuel_mass, pass, fail);

    % ---- Test 5: Better L/D → longer range (monotonicity) -------------
    state2 = state;
    state2.CL_CD = 20.0;   % better aerodynamics
    state2 = breguet(state2, dv, ac, atm, mis);
    [pass, fail] = check('Higher L/D → longer range', ...
        state2.range > state.range, pass, fail);

    % ---- Test 6: Lower TSFC → longer range (monotonicity) -------------
    state3 = state;
    state3.TSFC = 1.3e-5;   % more efficient engine
    state3 = breguet(state3, dv, ac, atm, mis);
    [pass, fail] = check('Lower TSFC → longer range', ...
        state3.range > state.range, pass, fail);

    % ---- Test 7: More fuel loaded → longer range -----------------------
    % Raise ac.fuel_mass directly; MTOW held fixed so climb penalty unchanged.
    ac2 = ac; ac2.fuel_mass = ac.fuel_mass * 1.10;   % 10% more fuel
    state4 = state;
    state4 = breguet(state4, dv, ac2, atm, mis);
    [pass, fail] = check('More fuel loaded (fuel_mass +10%) → longer range', ...
        state4.range > state.range, pass, fail);

    % ---- Test 8: Breguet formula sanity check (manual) -----------------
    % R = (V/(g*TSFC)) * L/D * ln(W_start/W_end) * ff_descent
    ff_to = 0.970; ff_climb = 0.985; ff_descent = 0.990;
    reserve_f = mis.reserve_f;
    W_start = ac.MTOW * ff_to * ff_climb;
    W_fuel_av = ac.fuel_mass * (1 - reserve_f) - ac.MTOW * (1 - ff_to * ff_climb);
    W_fuel_av = max(W_fuel_av, 1);
    W_end   = max(W_start - W_fuel_av, 0.5 * W_start);
    range_manual = (dv.V / (atm.g * state.TSFC)) * state.CL_CD * log(W_start / W_end) * ff_descent;
    [pass, fail] = check('Breguet output matches manual formula (within 0.1%)', ...
        abs(state.range - range_manual) / range_manual < 0.001, pass, fail);

    summary(pass, fail);
end

% =========================================================================
function [pass, fail] = check(name, condition, pass, fail)
    if condition
        fprintf('  [PASS] %s\n', name);
        pass = pass + 1;
    else
        fprintf('  [FAIL] %s\n', name);
        fail = fail + 1;
    end
end

function summary(pass, fail)
    fprintf('----------------------------------------\n');
    fprintf('  Results: %d passed, %d failed\n', pass, fail);
    if fail == 0
        fprintf('  ALL TESTS PASSED\n');
    else
        fprintf('  *** %d TEST(S) FAILED ***\n', fail);
    end
    fprintf('========================================\n\n');
end