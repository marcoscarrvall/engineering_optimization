clear all; clear global; clc; close all;

% --- Design Variables: [V, BPR, PR_fan, PR_LPC, PR_HPC] ---
x0 = [230, 8.0, 1.55, 1.55, 22.0];
lb = [200, 5.0, 1.1,  1.1,  8.0 ];
ub = [235, 14.0, 1.7, 2.5,  25.0];

data = A320data();

x0_norm  = normalize_vars(x0, lb, ub);
lb_norm  = zeros(size(x0_norm));
ub_norm  = ones(size(x0_norm));

mda_options.tol      = 1e-3;
mda_options.max_iter = 100;
mda_options.verbose  = false;

global optHistory
optHistory.iter   = [];
optHistory.fval   = [];
optHistory.x      = [];
optHistory.constr = [];
optHistory.count  = 0;

A = []; b = []; Aeq = []; beq = [];

optimizer_options = optimoptions('fmincon', ...
    'Algorithm',               'interior-point', ...
    'Display',                 'iter-detailed', ...
    'FiniteDifferenceType',    'central', ...
    'FiniteDifferenceStepSize', 5e-3, ...
    'OptimalityTolerance',     1e-6, ...
    'ConstraintTolerance',     1e-5, ...
    'StepTolerance',           1e-6, ...
    'OutputFcn',               @(x, optimValues, state) fmincon_history(x, optimValues, state, lb, ub, data, mda_options), ...
    'UseParallel',             false);

tic;
[x_opt_norm, f_opt, exitflag, output] = fmincon( ...
    @(x) optim(x, lb, ub, data, mda_options), ...
    x0_norm, A, b, Aeq, beq, lb_norm, ub_norm, ...
    @(x) constraints(x, lb, ub, data, mda_options), ...
    optimizer_options);
elapsed_time = toc;

x_opt = denormalize_vars(x_opt_norm, lb, ub);

var_names = {'V (m/s)', 'BPR', 'PR_{fan}', 'PR_{LPC}', 'PR_{HPC}'};

fprintf('\nOptimal Design Variables:\n');
for i = 1:length(x_opt)
    fprintf('  %s: %.3f\n', var_names{i}, x_opt(i));
end
fprintf('Optimal Objective: %.6f\n', f_opt);

% --- Convergence Plot ---
figure('Color', 'w', 'Name', 'Optimization Convergence');
plot(optHistory.iter, optHistory.fval, '-bo', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlim([0, max(optHistory.iter)]);
grid on;
ylabel('Objective f(x)');

% --- Design Variable Evolution ---
figure('Color', 'w', 'Name', 'Design Variable Evolution');
for v = 1:length(x0)
    subplot(length(x0), 1, v);
    hold on;

    yline(lb(v),   '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.2, 'Label', 'lb', 'LabelVerticalAlignment', 'bottom');
    yline(ub(v),   '--', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.2, 'Label', 'ub', 'LabelVerticalAlignment', 'top');
    yline(x0(v),   ':',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'Label', 'x_0', 'LabelVerticalAlignment', 'middle');
    yline(x_opt(v),'-',  'Color', [0.1 0.7 0.3], 'LineWidth', 1.5, 'Label', 'x_{opt}', 'LabelVerticalAlignment', 'middle');

    plot(optHistory.iter, optHistory.x(:, v), '-o', ...
        'LineWidth', 1.8, 'MarkerSize', 5, ...
        'Color', [0.1 0.45 0.85], 'MarkerFaceColor', [0.1 0.45 0.85]);

    y_range = optHistory.x(:, v);
    margin  = max((max(y_range) - min(y_range)) * 0.1, abs(min(y_range)) * 0.05 + 1e-6);
    ylim([min(y_range) - margin, max(y_range) + margin]);

    ylabel(var_names{v});
    grid on;
    if v == length(x0)
        xlabel('Iteration');
    else
        set(gca, 'XTickLabel', []);
    end
end

% --- Constraint History ---
constraint_names = {'Clearance', 'TIT Limit', 'Tip Mach'};

figure('Color', 'w', 'Name', 'Constraint History');
if ~isempty(optHistory.constr) && ~isempty(optHistory.iter)
    n      = min(length(optHistory.iter), size(optHistory.constr, 1));
    h      = plot(optHistory.iter(1:n), optHistory.constr(1:n, :), '-s', 'LineWidth', 1.2);
    hold on;
    yline(0, 'r--', 'LineWidth', 2);
    ylim([-1, 1]);
    grid on;
    ylabel('Value (c \leq 0)');
    xlabel('Evaluation Count');
    legend(h, constraint_names(1:size(optHistory.constr(1:n,:), 2)), 'Location', 'bestoutside');
end

% --- Flag MDA Failures (fval == 1) ---
tol_flag  = 1e-4;
idx_fail  = find(abs(optHistory.fval - 1) < tol_flag);

if ~isempty(idx_fail)
    fprintf('\nMDA failures detected (%d points with fval ≈ 1):\n', length(idx_fail));
    fprintf('%-4s | %-5s | %-5s | %-6s | %-6s | %-6s | Objective\n', ...
        'Iter', 'V', 'BPR', 'PR_fan', 'PR_LPC', 'PR_HPC');
    fprintf('%s\n', repmat('-', 1, 58));
    for i = 1:length(idx_fail)
        xf = optHistory.x(idx_fail(i), :);
        fprintf('%4d | %5.1f | %5.1f | %6.3f | %6.3f | %6.1f | %9.4f\n', ...
            idx_fail(i), xf(1), xf(2), xf(3), xf(4), xf(5), optHistory.fval(idx_fail(i)));
    end
else
    fprintf('\nNo MDA failures detected.\n');
end

fprintf('\nCompleted in %.2f s | Iterations: %d | Evaluations: %d\n', ...
    elapsed_time, output.iterations, output.funcCount);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function x_norm = normalize_vars(x, lb, ub)
    x_norm = (x - lb) ./ (ub - lb);
end

function x_real = denormalize_vars(x_norm, lb, ub)
    x_real = x_norm .* (ub - lb) + lb;
end

function stop = fmincon_history(x_norm, optimValues, state, lb, ub, data, mda_options)
    global optHistory
    stop = false;

    if strcmp(state, 'iter')
        optHistory.count = optHistory.count + 1;
        idx = optHistory.count;

        optHistory.iter(idx, 1) = optimValues.iteration;
        optHistory.fval(idx, 1) = optimValues.fval;
        optHistory.x(idx, :)    = x_norm(:)' .* (ub - lb) + lb;

        [g, ~] = constraints(x_norm, lb, ub, data, mda_options);
        optHistory.constr(idx, :) = g(:)';
    end
end


% =========================================================================
% ROBUSTNESS TEST: MULTIPLE STARTING POINTS
% =========================================================================
starting_points = [
    230, 8.0,  1.55, 1.55, 22.0;
    200, 5.0,  1.1,  1.1,  10.0;
    240, 14.0, 1.7,  2.0,  25.0;
    210, 6.5,  1.2,  1.3,  12.0;
    225, 11.0, 1.6,  1.5,  20.0;
    240, 14.0, 1.35, 1.8,  15.5;
];

n_starts    = size(starting_points, 1);
start_labels = {'Original', 'Lower Bounds', 'Upper Bounds', 'Feasible Low', 'Feasible Mid-High', 'Infeasible Candidate'};

results(n_starts) = struct('x0', [], 'x_opt', [], 'f_opt', [], 'g_opt', [], ...
                           'converged', false, 'feasible', false, 'iters', 0);

fprintf('\n--- Robustness Test (%d starting points) ---\n', n_starts);

for s = 1:n_starts
    fprintf('\n[Run %d/%d] %s\n', s, n_starts, start_labels{s});
    fprintf('  x0 = [%.1f, %.1f, %.2f, %.2f, %.1f]\n', starting_points(s,:));

    x_s = normalize_vars(starting_points(s,:), lb, ub);

    [x_opt_s, f_opt_s, exitflag_s, output_s] = fmincon( ...
        @(x) optim(x, lb, ub, data, mda_options), ...
        x_s, A, b, Aeq, beq, lb_norm, ub_norm, ...
        @(x) constraints(x, lb, ub, data, mda_options), ...
        optimizer_options);

    x_opt_phys = denormalize_vars(x_opt_s, lb, ub);   % <-- fixed: use x_opt_s, not x_s

    [g_opt, ~] = constraints(x_opt_s, lb, ub, data, mda_options);

    results(s).x0        = starting_points(s,:);
    results(s).x_opt     = x_opt_phys;
    results(s).f_opt     = f_opt_s;                    % <-- fixed: use fmincon output directly
    results(s).g_opt     = g_opt;
    results(s).converged = exitflag_s > 0;             % <-- fixed: use exitflag from fmincon
    results(s).feasible  = max(g_opt) <= 1e-4;
    results(s).iters     = output_s.iterations;        % <-- fixed: use output struct from fmincon
end

% =========================================================================
% RESULTS SUMMARY
% =========================================================================
fprintf('\n\n%s\n', repmat('=', 1, 95));
fprintf('ROBUSTNESS TEST SUMMARY\n');
fprintf('%s\n', repmat('=', 1, 95));
fprintf('%-22s  %-7s  %-7s  %-7s  %-7s  %-7s  %-10s  %-10s  %-5s  %-9s\n', ...
    'Start', 'V', 'BPR', 'PRfan', 'PRlpc', 'PRhpc', 'f_opt', 'max_g', 'Iters', 'Status');
fprintf('%s\n', repmat('-', 1, 95));

feasible_mask = [results.feasible];
f_opts        = [results.f_opt];
f_ref         = min(f_opts(feasible_mask));

for s = 1:n_starts
    r = results(s);
    if ~r.feasible
        status = 'INFEASIBLE';
    elseif ~r.converged
        status = 'NO CONV';
    elseif abs(r.f_opt - f_ref) < 1e-2
        status = 'GLOBAL';
    else
        status = 'LOCAL';
    end

    fprintf('%-22s  %-7.2f  %-7.2f  %-7.3f  %-7.3f  %-7.2f  %-10.6f  %-10.4f  %-5d  %-9s\n', ...
        start_labels{s}, r.x_opt(1), r.x_opt(2), r.x_opt(3), r.x_opt(4), r.x_opt(5), ...
        r.f_opt, max(r.g_opt), r.iters, status);
end

fprintf('\nBest feasible f_opt = %.6f\n', f_ref);
n_global = sum(arrayfun(@(r) r.feasible && r.converged && abs(r.f_opt - f_ref) < 1e-2, results));
fprintf('Runs at global minimum: %d / %d\n', n_global, n_starts);

if n_global == sum(feasible_mask & [results.converged])
    fprintf('Conclusion: Problem appears UNIMODAL.\n');
else
    fprintf('Conclusion: Problem appears MULTIMODAL — multiple local optima detected.\n');
end

% =========================================================================
% PLOT 1: DESIGN VARIABLE OUTCOMES
% =========================================================================
figure('Color', 'w', 'Name', 'Robustness: Final Design Variables');
x_matrix = reshape([results.x_opt], length(x0), n_starts)';
colors   = lines(n_starts);

for v = 1:length(x0)
    subplot(length(x0), 1, v); hold on;

    yline(lb(v),     '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.0);
    yline(ub(v),     '--', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.0);
    yline(x_opt(v),  '-',  'Color', [0.1 0.7 0.3], 'LineWidth', 1.5);

    for s = 1:n_starts
        if     ~results(s).feasible,  marker = 'x';
        elseif ~results(s).converged, marker = 'd';
        else,                         marker = 'o';
        end
        plot(s, x_matrix(s, v), marker, 'MarkerSize', 9, ...
            'Color', colors(s,:), 'MarkerFaceColor', colors(s,:), 'LineWidth', 1.5);
    end

    ylabel(var_names{v}); grid on;
    xlim([0.5, n_starts + 0.5]);
    set(gca, 'XTick', 1:n_starts);
    if v == length(x0)
        set(gca, 'XTickLabel', start_labels, 'XTickLabelRotation', 15);
    else
        set(gca, 'XTickLabel', []);
    end
end

% =========================================================================
% PLOT 2: OBJECTIVE VALUES
% =========================================================================
figure('Color', 'w', 'Name', 'Robustness: Objective Values');

bar_colors = zeros(n_starts, 3);
for s = 1:n_starts
    if ~results(s).feasible || ~results(s).converged
        bar_colors(s,:) = [0.8 0.2 0.2];           % red   — failed
    elseif abs(f_opts(s) - f_ref) < 1e-2
        bar_colors(s,:) = [0.1 0.7 0.3];           % green — global
    else
        bar_colors(s,:) = [1.0 0.6 0.0];           % amber — local
    end
end

b = bar(f_opts, 'FaceColor', 'flat');
b.CData = bar_colors;
yline(f_ref, 'k--', 'LineWidth', 1.5, 'Label', 'Global min', 'LabelVerticalAlignment', 'bottom');
set(gca, 'XTick', 1:n_starts, 'XTickLabel', start_labels, 'XTickLabelRotation', 15);
ylabel('f_{opt}');
grid on;