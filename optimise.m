function optimise()
%% OPTIMISE  —  MDO optimiser: minimise cruise range w.r.t. BPR
%
% Objective  :  minimise  state.range
% Design var :  BPR  in [BPR_lb, BPR_ub]
% Constraints (c <= 0):
%   g1 = TIT / TIT_max - 1          (turbine inlet temperature)
%   g2 = clearance_min/clearance - 1 (ground clearance)
%   g3 = M_tip / M_tip_max - 1      (fan-tip Mach / noise)
%
% thermo.m iterates FAR to match required thrust, so TIT varies with BPR.
%
% Solver: fmincon (interior-point) or fminsearch + penalty fallback.
%
% Usage:   optimise()
% Outputs: BPR_opt, state_opt, opt_hist  (written to base workspace)
% =========================================================================

run('setup.m');
[~, dv_ref, ~, con, ~, ~, ~, ~, ~] = build_reference();

BPR_lb = 4.0;
BPR_ub = 14.0;
BPR_0  = dv_ref.BPR;

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════╗\n');
fprintf('║   MDO OPTIMISER  —  Minimise Range  w.r.t.  BPR     ║\n');
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  BPR  in [%.1f, %.1f],  initial = %.2f               \n', BPR_lb, BPR_ub, BPR_0);
fprintf('║  Constraints : TIT, clearance, fan-tip Mach          \n');
fprintf('║  thermo      : FAR iterated to match thrust           \n');
fprintf('╚══════════════════════════════════════════════════════╝\n\n');

%% ---- Pre-check ---------------------------------------------------------
fprintf('Pre-check at BPR = %.2f:\n', BPR_0);
[f0, g0] = eval_mda(BPR_0, dv_ref, con);
fprintf('  Range      = %.1f km\n',    f0/1e3);
fprintf('  g_TIT      = %+.4f  %s\n', g0(1), feas(g0(1)));
fprintf('  g_clearance= %+.4f  %s\n', g0(2), feas(g0(2)));
fprintf('  g_M_tip    = %+.4f  %s\n', g0(3), feas(g0(3)));
fprintf('\n');

%% ---- Shared log --------------------------------------------------------
log    = {struct('BPR',[],'range',[],'g_TIT',[],'g_cl',[],'g_nse',[])};
n_eval = {0};

obj_fn = @(x) eval_obj(x, dv_ref, con, log, n_eval);
con_fn = @(x) eval_con(x, dv_ref, con);

%% ---- Solver ------------------------------------------------------------
use_fmincon = license('test','Optimization_Toolbox') && ~isempty(which('fmincon'));

if use_fmincon
    fprintf('Solver: fmincon (interior-point)\n\n');
    opts = optimoptions('fmincon', ...
        'Algorithm',            'interior-point', ...
        'Display',              'iter-detailed',  ...
        'FiniteDifferenceType', 'forward',        ...
        'FunctionTolerance',    1e-4,             ...
        'StepTolerance',        1e-4,             ...
        'MaxFunctionEvaluations', 300);
    [BPR_opt, range_opt, exitflag, output] = fmincon( ...
        obj_fn, BPR_0, [],[],[],[], BPR_lb, BPR_ub, con_fn, opts);
    fprintf('\n  Exit flag : %d  (%s)\n', exitflag, eflag(exitflag));
    fprintf('  Evals     : %d\n', output.funcCount);
else
    fprintf('No Optimisation Toolbox — fminsearch + penalty\n\n');
    mu     = 1e8;
    pen_fn = @(x) eval_penalty(x, mu, BPR_lb, BPR_ub, dv_ref, con, log, n_eval);
    opts_fm = optimset('Display','iter','TolX',1e-4,'TolFun',1e-4,'MaxFunEvals',500);
    [BPR_opt,~,exitflag,output] = fminsearch(pen_fn, max(BPR_lb,min(BPR_ub,BPR_0)), opts_fm);
    BPR_opt   = max(BPR_lb, min(BPR_ub, BPR_opt));
    range_opt = eval_obj(BPR_opt, dv_ref, con, log, n_eval);
    fprintf('\n  Exit flag : %d\n  Evals : %d\n', exitflag, output.funcCount);
end

%% ---- Result summary ----------------------------------------------------
fprintf('\n╔══════════════════════════════════════════════════════╗\n');
fprintf('║  OPTIMUM: BPR* = %.4f,  Range = %.1f km            \n', BPR_opt, range_opt/1e3);
fprintf('╚══════════════════════════════════════════════════════╝\n\n');

dv_opt     = dv_ref;
dv_opt.BPR = BPR_opt;
state_opt  = mda(dv_opt);

%% ---- Plots -------------------------------------------------------------
opt_hist = log{1};
n = numel(opt_hist.BPR);

if n > 0
    ev = 1:n;
    figure('Name','Optimiser History','NumberTitle','off','Position',[200 200 1000 600]);

    subplot(2,2,1);
    plot(ev, opt_hist.BPR,'b-o','LineWidth',1.5,'MarkerSize',5);
    xlabel('Evaluation'); ylabel('BPR [-]'); title('Design Variable'); grid on;

    subplot(2,2,2);
    plot(ev, opt_hist.range,'r-o','LineWidth',1.5,'MarkerSize',5);
    xlabel('Evaluation'); ylabel('Range [km]'); title('Objective (Range)'); grid on;

    subplot(2,2,3);
    h1 = plot(ev, opt_hist.g_TIT,'r-o','LineWidth',1.5,'MarkerSize',5); hold on;
    h2 = plot(ev, opt_hist.g_cl, 'b-s','LineWidth',1.5,'MarkerSize',5);
    h3 = plot(ev, opt_hist.g_nse,'g-^','LineWidth',1.5,'MarkerSize',5);
    hl = yline(0,'k--'); hl.HandleVisibility = 'off';
    legend([h1 h2 h3],'g_{TIT}','g_{clearance}','g_{M_{tip}}','Location','best');
    xlabel('Evaluation'); ylabel('g_i [-]'); title('Constraints (\leq0 feasible)'); grid on;

    subplot(2,2,4);
    scatter(opt_hist.BPR, opt_hist.range, 40, ev, 'filled');
    colorbar; xlabel('BPR [-]'); ylabel('Range [km]');
    title('Design-space exploration (colour = eval #)'); grid on;
    hl2 = xline(BPR_opt,'r--'); hl2.HandleVisibility = 'off';
    text(BPR_opt, min(opt_hist.range), sprintf(' BPR^*=%.2f',BPR_opt), ...
        'Color','r','FontSize',9,'VerticalAlignment','bottom');

    sgtitle('MDO Optimiser — Minimise Range w.r.t. BPR','FontSize',13,'FontWeight','bold');
else
    fprintf('No evaluations logged — skipping plots.\n');
end

assignin('base','BPR_opt',   BPR_opt);
assignin('base','state_opt', state_opt);
assignin('base','opt_hist',  opt_hist);
fprintf('\nResults saved to base workspace.\n');

end


%% =========================================================================
%%  LOCAL FUNCTIONS
%% =========================================================================

function [f, g] = eval_mda(BPR_val, dv_ref, con)
    dv = dv_ref; dv.BPR = BPR_val;
    state = mda_silent(dv);
    f    = state.range;
    g(1) = state.TIT         / con.TIT_max               - 1;
    g(2) = con.clearance_min / max(state.clearance, 1e-6) - 1;
    g(3) = state.M_tip       / con.M_tip_max              - 1;
end

function f = eval_obj(BPR_val, dv_ref, con, log, n_eval)
    [f, g] = eval_mda(BPR_val, dv_ref, con);
    n_eval{1} = n_eval{1} + 1;
    log{1}.BPR(end+1)   = BPR_val;
    log{1}.range(end+1) = f / 1e3;
    log{1}.g_TIT(end+1) = g(1);
    log{1}.g_cl(end+1)  = g(2);
    log{1}.g_nse(end+1) = g(3);
end

function [c, ceq] = eval_con(BPR_val, dv_ref, con)
    [~, g] = eval_mda(BPR_val, dv_ref, con);
    c = g(:); ceq = [];
end

function f_pen = eval_penalty(x, mu, lb, ub, dv_ref, con, log, n_eval)
    x = max(lb, min(ub, x));
    [f, g] = eval_mda(x, dv_ref, con);
    n_eval{1} = n_eval{1} + 1;
    log{1}.BPR(end+1)   = x;
    log{1}.range(end+1) = f / 1e3;
    log{1}.g_TIT(end+1) = g(1);
    log{1}.g_cl(end+1)  = g(2);
    log{1}.g_nse(end+1) = g(3);
    f_pen = f + mu*sum(max(0,g).^2);
end

function s = feas(g_val)
    if g_val <= 0; s = '(OK)'; else; s = '*** VIOLATED ***'; end
end

function s = eflag(flag)
    switch flag
        case  1;   s = 'optimality conditions met';
        case  2;   s = 'step size below TolX';
        case  0;   s = 'max iterations reached';
        case -1;   s = 'stopped by output function';
        case -2;   s = 'no feasible point found';
        otherwise; s = 'unknown';
    end
end