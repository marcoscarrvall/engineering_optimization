%% MDA_LOOP.m
% Multidisciplinary Analysis (MDA) Loop
%
% Couples WATE → AERO → THERMO → BREGUET in a fixed-point iteration until
% the required thrust (= total drag from AERO) converges.
%
% USAGE:
%   mda_loop          — runs with reference design vector (BPR=5.9, etc.)
%   mda_loop(dv)      — runs with a custom design vector dv
%
% EXAMPLE — sweep BPR:
%   [~, dv_ref, ~, ~, ~, ~, ~, ~, ~] = build_reference();
%   dv = dv_ref; dv.BPR = 10;
%   state = mda_loop(dv);
%
% OUTPUTS:
%   state   — converged state struct (all discipline outputs)
%   hist    — iteration history struct
%
% =========================================================================

function [state, hist] = mda(dv_in)

%% ---- 0.  Load reference data ------------------------------------------
[atm, dv_ref, ac, con, eng, thermo_data, ~, mis, state] = build_reference();

% Override design vector if one was supplied by the caller
if nargin == 1
    dv = dv_in;
else
    dv = dv_ref;
end

%% ---- 1.  MDA settings -------------------------------------------------
max_iter  = 50;
tol       = 1e-4;
omega     = 0.6;

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════╗\n');
fprintf('║   MDA LOOP  —  A320 / CFM56-5B                      ║\n');
fprintf('║   Coupling: WATE → AERO → THERMO → BREGUET          ║\n');
fprintf('╚══════════════════════════════════════════════════════╝\n');
fprintf('  BPR        : %.2f\n',   dv.BPR);
fprintf('  PR_Fan     : %.2f\n',   dv.PR_Fan);
fprintf('  PR_LPC     : %.2f\n',   dv.PR_LPC);
fprintf('  PR_HPC     : %.2f\n',   dv.PR_HPC);
fprintf('  V (m/s)    : %.1f\n',   dv.V);
fprintf('  Tolerance  : %.2e\n',   tol);
fprintf('  Relaxation : %.2f\n\n', omega);

% Iteration history arrays
hist.D_total  = zeros(1, max_iter);
hist.TSFC     = zeros(1, max_iter);
hist.mdot     = zeros(1, max_iter);
hist.W_engine = zeros(1, max_iter);
hist.MTOW     = zeros(1, max_iter);
hist.CL_CD    = zeros(1, max_iter);
hist.range    = zeros(1, max_iter);
hist.TIT      = zeros(1, max_iter);
hist.D_fan    = zeros(1, max_iter);

converged = false;
n_iter    = 0;

%% ---- 2.  Fixed-point iteration ----------------------------------------
for iter = 1:max_iter
    n_iter = iter;
    D_total_old = state.D_total;

    fprintf('--- Iteration %d -------------------------------------------\n', iter);

    state = wate(dv, atm, eng, state);
    state = aero(state, ac, atm, dv, eng);
    state = thermo(state, dv, atm, thermo_data);
    state = breguet(state, dv, ac, atm, mis);

    % Under-relaxation
    D_total_new   = state.D_total;
    state.D_total = omega * D_total_new + (1 - omega) * D_total_old;

    % Log
    hist.D_total(iter)  = state.D_total;
    hist.TSFC(iter)     = state.TSFC;
    hist.mdot(iter)     = state.mdot;
    hist.W_engine(iter) = state.W_engine;
    hist.MTOW(iter)     = state.MTOW;
    hist.CL_CD(iter)    = state.CL_CD;
    hist.range(iter)    = state.range;
    hist.TIT(iter)      = state.TIT;
    hist.D_fan(iter)    = state.D_fan;

    % Convergence check
    rel_change = abs(D_total_new - D_total_old) / max(abs(D_total_old), 1);
    fprintf('  D_total  (relaxed) = %10.2f N\n', state.D_total);
    fprintf('  Relative change    = %10.2e  (tol = %.2e)\n', rel_change, tol);

    if rel_change < tol && iter > 1
        converged = true;
        fprintf('\n  *** CONVERGED at iteration %d ***\n', iter);
        break;
    end
end

%% ---- 3.  Trim history arrays -----------------------------------------
hist.D_total  = hist.D_total(1:n_iter);
hist.TSFC     = hist.TSFC(1:n_iter);
hist.mdot     = hist.mdot(1:n_iter);
hist.W_engine = hist.W_engine(1:n_iter);
hist.MTOW     = hist.MTOW(1:n_iter);
hist.CL_CD    = hist.CL_CD(1:n_iter);
hist.range    = hist.range(1:n_iter);
hist.TIT      = hist.TIT(1:n_iter);
hist.D_fan    = hist.D_fan(1:n_iter);

%% ---- 4.  Converged solution summary ----------------------------------
fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════╗\n');
if converged
    fprintf('║   CONVERGED SOLUTION  (%d iterations)               \n', n_iter);
else
    fprintf('║   *** DID NOT CONVERGE in %d iterations ***         \n', n_iter);
end
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  Engine & Geometry                                   ║\n');
fprintf('║    D_fan         = %7.4f m                          \n', state.D_fan);
fprintf('║    D_nacelle     = %7.4f m                          \n', state.D_nacelle);
fprintf('║    L_engine      = %7.4f m                          \n', state.L_engine);
fprintf('║    W_engine      = %7.1f kg  (per engine)           \n', state.W_engine);
fprintf('║    mdot          = %7.2f kg/s (both engines)        \n', state.mdot);
fprintf('║    N_stages      = %7d                              \n', state.N_stages);
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  Thermodynamics                                      ║\n');
fprintf('║    TIT           = %7.1f K                          \n', state.TIT);
fprintf('║    TSFC          = %.4e kg/N/s                    \n', state.TSFC);
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  Aerodynamics & Weights                              ║\n');
fprintf('║    MTOW          = %7.1f kg                         \n', state.MTOW);
fprintf('║    CL            = %7.5f                            \n', state.CL);
fprintf('║    L/D           = %7.4f                            \n', state.CL_CD);
fprintf('║    D_total       = %7.1f N                          \n', state.D_total);
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  Mission                                             ║\n');
fprintf('║    Range         = %7.1f km  (%6.1f nm)            \n', state.range/1e3, state.range/1852);
fprintf('║    Design range  = %7.1f km                         \n', mis.R_design/1e3);
fprintf('╚══════════════════════════════════════════════════════╝\n\n');

%% ---- 5.  Constraint evaluation ---------------------------------------
fprintf('--- Constraint Check (converged state) ---\n');
[c_TIT, v_TIT, m_TIT]       = tit(state, con);
[c_cl,  v_cl,  m_cl]        = clearance(state, eng, con);
[c_nse, v_nse, m_nse, state] = noise(state, dv, atm, thermo_data, con);

fprintf('\n  Summary:\n');
fprintf('    TIT       : c = %+.6f,  margin = %+.6f  %s\n', c_TIT, m_TIT, status_str(v_TIT));
fprintf('    Clearance : c = %+.6f,  margin = %+.6f  %s\n', c_cl,  m_cl,  status_str(v_cl));
fprintf('    Noise     : c = %+.6f,  margin = %+.6f  %s\n', c_nse, m_nse, status_str(v_nse));

%% ---- 6.  Convergence plot --------------------------------------------
iters = 1:n_iter;

figure('Name', sprintf('MDA  BPR=%.1f', dv.BPR), 'NumberTitle', 'off', 'Position', [100 100 1100 700]);

subplot(3,3,1); plot(iters, hist.D_total/1e3,'b-o','LineWidth',1.5,'MarkerSize',5);
xlabel('Iteration'); ylabel('D_{total} [kN]'); title('Required Thrust'); grid on;

subplot(3,3,2); plot(iters, hist.mdot,'r-o','LineWidth',1.5,'MarkerSize',5);
xlabel('Iteration'); ylabel('\dot{m} [kg/s]'); title('Total Mass Flow'); grid on;

subplot(3,3,3); plot(iters, hist.W_engine,'g-o','LineWidth',1.5,'MarkerSize',5);
xlabel('Iteration'); ylabel('W_{engine} [kg]'); title('Engine Weight (per eng)'); grid on;

subplot(3,3,4); plot(iters, hist.MTOW,'m-o','LineWidth',1.5,'MarkerSize',5);
xlabel('Iteration'); ylabel('MTOW [kg]'); title('MTOW'); grid on;

subplot(3,3,5); plot(iters, hist.CL_CD,'k-o','LineWidth',1.5,'MarkerSize',5);
xlabel('Iteration'); ylabel('L/D [-]'); title('Lift-to-Drag Ratio'); grid on;

subplot(3,3,6); plot(iters, hist.TSFC*1e6,'c-o','LineWidth',1.5,'MarkerSize',5);
xlabel('Iteration'); ylabel('TSFC [mg/N/s]'); title('TSFC'); grid on;

subplot(3,3,7); plot(iters, hist.TIT,'r-o','LineWidth',1.5,'MarkerSize',5);
xlabel('Iteration'); ylabel('TIT [K]'); title('Turbine Inlet Temperature'); grid on;

subplot(3,3,8); plot(iters, hist.D_fan,'b-o','LineWidth',1.5,'MarkerSize',5);
xlabel('Iteration'); ylabel('D_{fan} [m]'); title('Fan Diameter'); grid on;

subplot(3,3,9); plot(iters, hist.range/1e3,'g-o','LineWidth',1.5,'MarkerSize',5);
xlabel('Iteration'); ylabel('Range [km]'); title('Cruise Range'); grid on;
yline(mis.R_design/1e3,'r--','Design Range','LabelHorizontalAlignment','left');

sgtitle(sprintf('MDA Convergence  —  BPR=%.1f  (%s,  %d iters)', ...
    dv.BPR, convergence_label(converged), n_iter), 'FontSize', 13, 'FontWeight', 'bold');

end

%% ---- Helper functions ------------------------------------------------
function s = status_str(violated)
    if violated; s = '<-- VIOLATED'; else; s = 'OK'; end
end

function s = convergence_label(converged)
    if converged; s = 'CONVERGED'; else; s = 'NOT CONVERGED'; end
end