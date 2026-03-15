%% TEST_THERMO.m
% Unit tests for THERMO.m
%
% Tests that the 0-D cycle analysis produces physically consistent outputs
% for the A320/CFM56-5B reference case.
% =========================================================================

function test_thermo()
    fprintf('========================================\n');
    fprintf(' TEST SUITE: THERMO.m\n');
    fprintf('========================================\n');
    pass = 0; fail = 0;

    % ---- Build reference structs ----------------------------------------
    [atm, dv, ~, ~, ~, thermo_data, ~, ~, state] = build_reference();

    % State requires D_total (total drag = required thrust)
    state.D_total = 120000;   % [N] ~60 kN per engine, reasonable cruise thrust

    % ---- Run module -------------------------------------------------------
    state = thermo(state, dv, atm, thermo_data);

    % ---- Test 1: TIT in physically plausible range ----------------------
    [pass, fail] = check('TIT in range [1200 K, 2000 K]', ...
        state.TIT >= 1200 && state.TIT <= 2000, pass, fail);

    % ---- Test 2: TIT matches FAR-derived expectation --------------------
    % T04 = T_HPC_exit + FAR*eta_cc*LHV / Cp_h
    % We can't easily recompute T03 here, but TIT should be >> 1000 K
    [pass, fail] = check('TIT > 1400 K (sanity for cruise FAR=0.027)', ...
        state.TIT > 1400, pass, fail);

    % ---- Test 3: TSFC in realistic range --------------------------------
    % Typical cruise TSFC for CFM56-5B ~1.5e-5 to 2.0e-5 kg/N/s
    [pass, fail] = check('TSFC in range [1e-5, 3e-5] kg/N/s', ...
        state.TSFC >= 1e-5 && state.TSFC <= 3e-5, pass, fail);

    % ---- Test 4: mdot positive and reasonable ---------------------------
    % CFM56-5B total flow ~300-400 kg/s for both engines at cruise power
    [pass, fail] = check('mdot (both engines) > 0', ...
        state.mdot > 0, pass, fail);
    [pass, fail] = check('mdot (both engines) in range [50, 800] kg/s', ...
        state.mdot >= 50 && state.mdot <= 800, pass, fail);

    % ---- Test 5: Specific thrust positive -------------------------------
    Fsp = state.D_total / state.mdot;
    [pass, fail] = check('Specific thrust Fsp > 0', Fsp > 0, pass, fail);

    % ---- Test 6: TSFC-Fsp consistency  (TSFC = FAR/((1+BPR)*Fsp)) -----
    FAR = thermo_data.FAR; BPR = dv.BPR;
    TSFC_check = (FAR / (1 + BPR)) / Fsp;
    [pass, fail] = check('TSFC consistent with Fsp (within 1%)', ...
        abs(state.TSFC - TSFC_check) / TSFC_check < 0.01, pass, fail);

    % ---- Test 7: Higher thrust → higher mdot ----------------------------
    state2 = state; state2.D_total = 200000;
    state2 = thermo(state2, dv, atm, thermo_data);
    [pass, fail] = check('Higher thrust demand → higher mdot', ...
        state2.mdot > state.mdot, pass, fail);

    % ---- Test 8: Higher BPR → lower TSFC (efficiency benefit) ----------
    dv2 = dv; dv2.BPR = 10;
    state3 = state; state3.D_total = 120000;
    state3 = thermo(state3, dv2, atm, thermo_data);
    [pass, fail] = check('Higher BPR → lower TSFC', ...
        state3.TSFC < state.TSFC, pass, fail);

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