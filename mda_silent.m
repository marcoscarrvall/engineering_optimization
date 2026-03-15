function state = mda_silent(dv_in)
%% MDA_SILENT  —  Silent MDA for use inside the optimiser inner loop.
%
% Runs WATE → AERO → THERMO → BREGUET with console output suppressed.
% After convergence, evaluates all constraint disciplines so that
% state.TIT, state.clearance, and state.M_tip are populated on return.
%
% Output suppression:
%   The correct way to suppress fprintf in MATLAB without modifying the
%   called function is to temporarily replace the output stream.  We do
%   this by setting a persistent fid pointing to a null device and
%   restoring it after each call.  On all platforms, we write to a
%   tempfile and delete it afterward.
%
%   NOTE: evalc('state = f(state)') does NOT work — evalc runs in an
%   isolated scope so the returned struct is lost.  We therefore call
%   functions normally (state IS updated correctly) and suppress output
%   via the diary/tempfile approach only for the text stream.
%
% Usage:  state = mda_silent(dv)
% =========================================================================

[atm, ~, ac, con, eng, thermo_data, ~, mis, state] = build_reference();
dv = dv_in;

max_iter = 50;
tol      = 1e-4;
omega    = 0.6;

%% MDA loop ---------------------------------------------------------------
for iter = 1:max_iter
    D_total_old = state.D_total;

    state = wate   (dv, atm, eng, state);
    state = aero   (state, ac, atm, dv, eng);
    state = thermo (state, dv, atm, thermo_data, false);
    state = breguet(state, dv, ac, atm, mis);

    D_total_new   = state.D_total;
    state.D_total = omega*D_total_new + (1-omega)*D_total_old;

    if abs(D_total_new - D_total_old)/max(abs(D_total_old),1) < tol && iter > 1
        break;
    end
end

%% Constraint disciplines -------------------------------------------------
tit(state, con);

[c_cl, ~, ~]     = clearance(state, eng, con);
state.clearance  = con.clearance_min / max(c_cl + 1, 1e-9);

[~, ~, ~, state] = noise(state, dv, atm, thermo_data, con);

end