%% TEST_AERO.m
% Unit tests for AERO.m
%
% Tests that the aerodynamics module correctly updates MTOW, computes CL,
% drag increments, and lift-to-drag ratio.
% =========================================================================

function test_aero()
    fprintf('========================================\n');
    fprintf(' TEST SUITE: AERO.m\n');
    fprintf('========================================\n');
    pass = 0; fail = 0;

    % ---- Build reference structs ----------------------------------------
    [atm, dv, ac, ~, eng, ~, ~, ~, state] = build_reference();

    % Provide engine geometry (normally from WATE)
    state.W_engine  = eng.W_engine_ref;   % [kg]  at reference → no MTOW delta
    state.D_fan     = eng.D_fan_ref;      % [m]
    state.D_nacelle = eng.D_fan_ref * 1.10;
    state.L_engine  = eng.L_engine_ref;   % [m]

    % ---- Run module -------------------------------------------------------
    state = aero(state, ac, atm, dv, eng);

    % ---- Test 1: MTOW is positive and bounded ---------------------------
    [pass, fail] = check('MTOW > 0', state.MTOW > 0, pass, fail);
    [pass, fail] = check('MTOW >= 0.80 * ac.MTOW (lower bound enforced)', ...
        state.MTOW >= 0.80 * ac.MTOW, pass, fail);

    % ---- Test 2: At reference engine weight, MTOW = OEW+payload+fuel ---
    MTOW_expected = ac.OEW + ac.payload + ac.fuel_mass;
    [pass, fail] = check('Reference engine weight → MTOW = OEW+payload+fuel', ...
        abs(state.MTOW - MTOW_expected) < 1e-6, pass, fail);

    % ---- Test 3: CL in physically realistic range -----------------------
    [pass, fail] = check('CL > 0', state.CL > 0, pass, fail);
    [pass, fail] = check('CL in range [0.1, 1.5]', ...
        state.CL >= 0.1 && state.CL <= 1.5, pass, fail);
    % Close to reference cruise CL ~ 0.56
    [pass, fail] = check('CL within 20% of reference CL (0.56)', ...
        abs(state.CL - ac.CL_cr) / ac.CL_cr < 0.20, pass, fail);

    % ---- Test 4: L/D positive and in plausible range -------------------
    [pass, fail] = check('L/D > 0', state.CL_CD > 0, pass, fail);
    [pass, fail] = check('L/D in range [5, 30]', ...
        state.CL_CD >= 5 && state.CL_CD <= 30, pass, fail);

    % ---- Test 5: Total drag positive and sensible ----------------------
    [pass, fail] = check('D_total > 0', state.D_total > 0, pass, fail);
    [pass, fail] = check('D_total in range [20 kN, 500 kN]', ...
        state.D_total >= 20e3 && state.D_total <= 500e3, pass, fail);

    % ---- Test 6: Heavier engine → higher MTOW --------------------------
    state2 = state;
    state2.W_engine = eng.W_engine_ref + 500;  % 500 kg heavier per engine
    state2 = aero(state2, ac, atm, dv, eng);
    [pass, fail] = check('Heavier engine → higher MTOW', ...
        state2.MTOW > state.MTOW, pass, fail);

    % ---- Test 7: Larger nacelle → more drag, lower L/D -----------------
    state3 = state;
    state3.D_fan     = eng.D_fan_ref * 1.30;   % 30% bigger fan
    state3.D_nacelle = state3.D_fan * 1.10;
    state3 = aero(state3, ac, atm, dv, eng);
    [pass, fail] = check('Larger nacelle → lower L/D', ...
        state3.CL_CD < state.CL_CD, pass, fail);

    % ---- Test 8: L/D = CL / CD_total consistency -----------------------
    q   = 0.5 * atm.rho_cr * dv.V^2;
    CL_check = (state.MTOW * atm.g) / (q * ac.S_ref);
    [pass, fail] = check('CL consistent with level flight (L=W, within 0.1%)', ...
        abs(state.CL - CL_check) / CL_check < 0.001, pass, fail);

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