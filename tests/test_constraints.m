%% TEST_CONSTRAINTS.m
% Unit tests for CLEARANCE.m, NOISE.m, and TIT.m
%
% All three constraints now return NORMALISED (dimensionless) values:
%   c  = (value - limit) / limit   <= 0
%   margin = -c  (positive = safe headroom, expressed as a fraction of limit)
% =========================================================================

function test_constraints()
    test_clearance();
    test_noise();
    test_tit();
end

% =========================================================================
%  CLEARANCE
% =========================================================================
function test_clearance()
    fprintf('========================================\n');
    fprintf(' TEST SUITE: CLEARANCE.m\n');
    fprintf('========================================\n');
    pass = 0; fail = 0;

    [~, ~, ~, con, eng, ~, ~, ~, state] = build_reference();

    % ---- Case A: Reference nacelle — should be satisfied ----------------
    state.D_nacelle = 1.735 * 1.10;   % D_nacelle_ref = 1.909 m

    [c_cl, violated, margin] = clearance(state, eng, con);

    cl_val = eng.wing_h_AGL - eng.pylon_h - state.D_nacelle/2;

    [pass, fail] = check('[Ref] physical clearance > 0', cl_val > 0, pass, fail);

    % c_cl = (clearance_min - clearance) / clearance_min
    c_expected = (con.clearance_min - cl_val) / con.clearance_min;
    [pass, fail] = check('[Ref] c_cl matches normalised formula', ...
        abs(c_cl - c_expected) < 1e-9, pass, fail);

    % margin = -c_cl
    [pass, fail] = check('[Ref] margin = -c_cl', ...
        abs(margin + c_cl) < 1e-9, pass, fail);

    [pass, fail] = check('[Ref] violated == (c_cl > 0)', ...
        violated == (c_cl > 0), pass, fail);

    % c_cl dimensionless: |c_cl| should be O(1) or smaller
    [pass, fail] = check('[Ref] c_cl is dimensionless O(1)', ...
        abs(c_cl) < 10, pass, fail);

    % ---- Case B: Oversized nacelle — VIOLATED ---------------------------
    state.D_nacelle = 3.0;
    [c_cl2, violated2, margin2] = clearance(state, eng, con);
    [pass, fail] = check('[Oversized] constraint VIOLATED', violated2 == true, pass, fail);
    [pass, fail] = check('[Oversized] c_cl > 0', c_cl2 > 0, pass, fail);
    [pass, fail] = check('[Oversized] margin < 0', margin2 < 0, pass, fail);

    % ---- Case C: Small nacelle — SATISFIED ------------------------------
    state.D_nacelle = 1.0;
    [c_cl3, violated3, margin3] = clearance(state, eng, con);
    [pass, fail] = check('[Small nac] constraint SATISFIED', violated3 == false, pass, fail);
    [pass, fail] = check('[Small nac] c_cl <= 0', c_cl3 <= 0, pass, fail);
    [pass, fail] = check('[Small nac] margin > 0', margin3 > 0, pass, fail);

    % ---- Case D: Monotonicity — larger nacelle → tighter (larger c_cl) --
    state.D_nacelle = 1.5;
    [c_cl_small, ~, ~] = clearance(state, eng, con);
    state.D_nacelle = 2.0;
    [c_cl_large, ~, ~] = clearance(state, eng, con);
    [pass, fail] = check('Larger D_nacelle → larger c_cl (tighter)', ...
        c_cl_large > c_cl_small, pass, fail);

    % ---- Case E: At exact limit — c_cl = 0 -----------------------------
    D_at_limit = 2 * (eng.wing_h_AGL - eng.pylon_h - con.clearance_min);
    state.D_nacelle = D_at_limit;
    [c_cl_lim, ~, margin_lim] = clearance(state, eng, con);
    [pass, fail] = check('[At limit] c_cl = 0', abs(c_cl_lim) < 1e-9, pass, fail);
    [pass, fail] = check('[At limit] margin = 0', abs(margin_lim) < 1e-9, pass, fail);

    summary(pass, fail);
end

% =========================================================================
%  NOISE
% =========================================================================
function test_noise()
    fprintf('========================================\n');
    fprintf(' TEST SUITE: NOISE.m\n');
    fprintf('========================================\n');
    pass = 0; fail = 0;

    [atm, dv, ~, con, eng, thermo_data, ~, ~, state] = build_reference();

    state.D_fan = eng.D_fan_ref;
    state.A_fan = pi/4 * eng.D_fan_ref^2 * (1 - 0.30^2);

    [c_noise, violated, margin, state] = noise(state, dv, atm, thermo_data, con);

    % ---- Test 1: M_tip positive and in range ----------------------------
    [pass, fail] = check('M_tip > 0', state.M_tip > 0, pass, fail);
    [pass, fail] = check('M_tip in range [0.8, 2.0]', ...
        state.M_tip >= 0.8 && state.M_tip <= 2.0, pass, fail);

    % ---- Test 2: Normalised formula c = (M_tip - M_tip_max)/M_tip_max --
    c_expected = (state.M_tip - con.M_tip_max) / con.M_tip_max;
    [pass, fail] = check('c_noise matches normalised formula', ...
        abs(c_noise - c_expected) < 1e-9, pass, fail);

    % ---- Test 3: margin = -c_noise -------------------------------------
    [pass, fail] = check('margin = -c_noise', ...
        abs(margin + c_noise) < 1e-9, pass, fail);

    % ---- Test 4: violated flag consistent ------------------------------
    [pass, fail] = check('violated flag consistent with c_noise', ...
        violated == (c_noise > 0), pass, fail);

    % ---- Test 5: c_noise dimensionless O(1) ----------------------------
    [pass, fail] = check('c_noise is dimensionless O(1)', ...
        abs(c_noise) < 5, pass, fail);

    % ---- Test 6: Higher fan PR → larger c_noise (tighter) --------------
    dv_hi = dv; dv_hi.PR_Fan = 2.5;
    [c_noise_hi, ~, ~, state_hi] = noise(state, dv_hi, atm, thermo_data, con);
    [pass, fail] = check('Higher fan PR → higher M_tip', ...
        state_hi.M_tip > state.M_tip, pass, fail);
    [pass, fail] = check('Higher fan PR → larger c_noise', ...
        c_noise_hi > c_noise, pass, fail);

    % ---- Test 7: Tight limit triggers violation ------------------------
    con_tight = con; con_tight.M_tip_max = 0.5;
    [~, viol_tight, ~, ~] = noise(state, dv, atm, thermo_data, con_tight);
    [pass, fail] = check('Tight limit (M_tip_max=0.5) triggers violation', ...
        viol_tight == true, pass, fail);

    % ---- Test 8: At exact limit — c_noise = 0 -------------------------
    con_exact = con; con_exact.M_tip_max = state.M_tip;
    [c_exact, ~, margin_exact, ~] = noise(state, dv, atm, thermo_data, con_exact);
    [pass, fail] = check('[At limit] c_noise = 0', abs(c_exact) < 1e-9, pass, fail);
    [pass, fail] = check('[At limit] margin = 0', abs(margin_exact) < 1e-9, pass, fail);

    summary(pass, fail);
end

% =========================================================================
%  TIT
% =========================================================================
function test_tit()
    fprintf('========================================\n');
    fprintf(' TEST SUITE: TIT.m\n');
    fprintf('========================================\n');
    pass = 0; fail = 0;

    [~, ~, ~, con, ~, ~, ~, ~, state] = build_reference();

    % ---- Case A: TIT well below limit — SATISFIED -----------------------
    state.TIT = 1480;
    [c_TIT, violated, margin] = tit(state, con);
    [pass, fail] = check('[1480 K] constraint SATISFIED', violated == false, pass, fail);
    [pass, fail] = check('[1480 K] c_TIT <= 0', c_TIT <= 0, pass, fail);
    [pass, fail] = check('[1480 K] margin > 0', margin > 0, pass, fail);

    % ---- Case B: TIT exactly at limit — borderline ----------------------
    state.TIT = con.TIT_max;
    [c_TIT2, violated2, margin2] = tit(state, con);
    [pass, fail] = check('[At limit] c_TIT = 0', abs(c_TIT2) < 1e-9, pass, fail);
    [pass, fail] = check('[At limit] violated = false', violated2 == false, pass, fail);
    [pass, fail] = check('[At limit] margin = 0', abs(margin2) < 1e-9, pass, fail);

    % ---- Case C: TIT above limit — VIOLATED -----------------------------
    state.TIT = con.TIT_max + 50;
    [c_TIT3, violated3, margin3] = tit(state, con);
    [pass, fail] = check('[+50 K] constraint VIOLATED', violated3 == true, pass, fail);
    c_expected3 = 50 / con.TIT_max;
    [pass, fail] = check('[+50 K] c_TIT = 50/TIT_max (normalised)', ...
        abs(c_TIT3 - c_expected3) < 1e-9, pass, fail);
    [pass, fail] = check('[+50 K] margin = -c_TIT', ...
        abs(margin3 + c_TIT3) < 1e-9, pass, fail);

    % ---- Case D: Normalised formula and margin consistency --------------
    state.TIT = 1600;
    [c_TIT4, ~, margin4] = tit(state, con);
    c_expected4 = (1600 - con.TIT_max) / con.TIT_max;
    [pass, fail] = check('c_TIT = (TIT - TIT_max)/TIT_max', ...
        abs(c_TIT4 - c_expected4) < 1e-9, pass, fail);
    [pass, fail] = check('margin = -c_TIT', ...
        abs(margin4 + c_TIT4) < 1e-9, pass, fail);

    % ---- Case E: c_TIT dimensionless O(1) for typical TIT range --------
    state.TIT = 1900;
    [c_TIT5, ~, ~] = tit(state, con);
    [pass, fail] = check('c_TIT dimensionless O(1) even at 1900 K', ...
        abs(c_TIT5) < 2, pass, fail);

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