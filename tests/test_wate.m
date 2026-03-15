%% TEST_WATE.m
% Unit tests for WATE.m
%
% Tests that the engine weight and geometry sizing produces physically
% consistent outputs for the A320/CFM56-5B reference case.
% =========================================================================

function test_wate()
    fprintf('========================================\n');
    fprintf(' TEST SUITE: WATE.m\n');
    fprintf('========================================\n');
    pass = 0; fail = 0;

    % ---- Build reference structs ----------------------------------------
    [atm, dv, ~, ~, eng, ~, ~, ~, state] = build_reference();

    state.mdot = 350;   % [kg/s] total for both engines (reference)

    % ---- Run module -------------------------------------------------------
    state = wate(dv, atm, eng, state);

    % ---- Test 1: Fan diameter close to CFM56-5B reference (1.735 m) -----
    [pass, fail] = check('D_fan > 0', state.D_fan > 0, pass, fail);
    [pass, fail] = check('D_fan in range [0.8, 3.0] m', ...
        state.D_fan >= 0.8 && state.D_fan <= 3.0, pass, fail);
    [pass, fail] = check('D_fan within 30% of CFM56-5B ref (1.735 m)', ...
        abs(state.D_fan - 1.735) / 1.735 < 0.30, pass, fail);

    % ---- Test 2: Nacelle = 1.10 × D_fan --------------------------------
    [pass, fail] = check('D_nacelle = 1.10 * D_fan', ...
        abs(state.D_nacelle - 1.10 * state.D_fan) < 1e-9, pass, fail);

    % ---- Test 3: Engine length positive and sensible --------------------
    [pass, fail] = check('L_engine > 0', state.L_engine > 0, pass, fail);
    [pass, fail] = check('L_engine in range [1.5, 8.0] m', ...
        state.L_engine >= 1.5 && state.L_engine <= 8.0, pass, fail);

    % ---- Test 4: Fan annulus area consistent with continuity ------------
    % A_fan should equal mdot_per_eng / (rho_ax * V_ax) from continuity
    mdot_per_eng = state.mdot / 2;
    M0   = dv.V / atm.a_cr;
    T02  = atm.T_cr * (1 + (atm.gamma-1)/2 * M0^2);
    P02  = atm.P_cr * (T02 / atm.T_cr)^(atm.gamma / (atm.gamma-1));
    M_ax = 0.55;
    T_ax = T02 / (1 + (atm.gamma-1)/2 * M_ax^2);
    P_ax = P02 / (1 + (atm.gamma-1)/2 * M_ax^2)^(atm.gamma / (atm.gamma-1));
    rho_ax = P_ax / (atm.R * T_ax);
    V_ax = M_ax * sqrt(atm.gamma * atm.R * T_ax);
    A_fan_expected = mdot_per_eng / (rho_ax * V_ax);
    [pass, fail] = check('A_fan consistent with continuity (within 1%)', ...
        abs(state.A_fan - A_fan_expected) / A_fan_expected < 0.01, pass, fail);

    % ---- Test 5: Fan area consistent with fan diameter ------------------
    hub_tip = 0.30;
    D_fan_from_area = sqrt(4 * state.A_fan / (pi * (1 - hub_tip^2)));
    [pass, fail] = check('D_fan consistent with A_fan (within 0.1%)', ...
        abs(state.D_fan - D_fan_from_area) / D_fan_from_area < 0.001, pass, fail);

    % ---- Test 6: Stage count reasonable for OPR -------------------------
    OPR = dv.PR_Fan * dv.PR_LPC * dv.PR_HPC;
    N_expected = ceil(log(OPR) / log(1.25));
    [pass, fail] = check('N_stages matches ceil(log(OPR)/log(1.25))', ...
        state.N_stages == N_expected, pass, fail);

    % ---- Test 7: Engine weight positive and sensible --------------------
    [pass, fail] = check('W_engine > 0', state.W_engine > 0, pass, fail);
    [pass, fail] = check('W_engine in range [1000, 8000] kg', ...
        state.W_engine >= 1000 && state.W_engine <= 8000, pass, fail);

    % ---- Test 8: Higher mdot → larger fan, higher weight ----------------
    state2 = state; state2.mdot = 700;
    state2 = wate(dv, atm, eng, state2);
    [pass, fail] = check('Higher mdot → larger D_fan', ...
        state2.D_fan > state.D_fan, pass, fail);
    [pass, fail] = check('Higher mdot → higher W_engine', ...
        state2.W_engine > state.W_engine, pass, fail);

    % ---- Test 9: Higher BPR → bigger fan (same thrust demand) ----------
    % Physics: higher BPR → lower specific thrust (Fsp) → more mdot needed
    % for the same thrust. D_fan depends on mdot, not BPR directly.
    % Simulate the extra mdot a high-BPR engine needs for identical thrust.
    % Rough scaling: mdot ∝ (1+BPR) at fixed core flow.
    dv2 = dv; dv2.BPR = 10;
    mdot_scaled = state.mdot * (1 + dv2.BPR) / (1 + dv.BPR);
    state3 = state; state3.mdot = mdot_scaled;
    state3 = wate(dv2, atm, eng, state3);
    [pass, fail] = check('Higher BPR + scaled mdot → larger D_fan', ...
        state3.D_fan > state.D_fan, pass, fail);

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