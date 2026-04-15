clear; clc; close all;

% =========================================================================
% PROBLEM SETUP
% =========================================================================
x0 = [230, 8.0];
lb = [200, 5.0];
ub = [240, 14.0];

var_names        = {'V (m/s)', 'BPR'};
constraint_names = {'Clearance', 'TIT Limit', 'Tip Mach'};

global optHistory
optHistory.iter   = [];
optHistory.fval   = [];
optHistory.x      = [];
optHistory.constr = [];

global history
history.iter   = [];
history.fval   = [];
history.x      = [];
history.constr = [];
history.count  = 0;

mda_options.tol      = 1e-3;
mda_options.max_iter = 100;
mda_options.verbose  = false;
data = A320data;

% =========================================================================
% SLP SETTINGS
% =========================================================================
max_slp_iter = 100;
move_limit   = 0.1;
M_penalty    = 1e6;
h_step       = 5e-3;
conv_tol     = 1e-4;

x_curr      = normalize_vars(x0, lb, ub);
multipliers = [];

fprintf('\n--- Starting SLP Optimization ---\n');

% =========================================================================
% MAIN LOOP
% =========================================================================
for k = 1:max_slp_iter

    f_curr = simplified_optim(x_curr, lb, ub, data, mda_options);
    [g_curr, ~] = simplified_constraints(x_curr, lb, ub, data, mda_options);

    history.count        = history.count + 1;
    history.iter(k)      = k;
    history.fval(k)      = f_curr;
    history.x(k, :)      = denormalize_vars(x_curr, lb, ub);
    history.constr(k, :) = g_curr';

    [grad_f, grad_g] = get_gradients(x_curr, lb, ub, data, mda_options, h_step);
    grad_f = grad_f(:);

    f_lin         = [grad_f; M_penalty];
    n_constraints = size(grad_g, 1);
    A_lin         = [grad_g, -ones(n_constraints, 1)];
    b_lin         = -g_curr;

    lb_dx  = max(-move_limit, 0 - x_curr);
    ub_dx  = min(move_limit,  1 - x_curr);
    lb_lin = [lb_dx(:); 0];
    ub_lin = [ub_dx(:); inf];

    options_lp = optimoptions('linprog', 'Display', 'none');
    [res, ~, exitflag, ~, lambda_struct] = linprog(f_lin, A_lin, b_lin, [], [], lb_lin, ub_lin, options_lp);

    if exitflag ~= 1
        fprintf('Iter %d: linprog did not converge. Stopping.\n', k); break;
    end

    dx          = res(1:2)';
    beta        = res(3);
    multipliers = lambda_struct.ineqlin;

    x_next = x_curr + dx;
    f_next = simplified_optim(x_next, lb, ub, data, mda_options);
    [g_next, ~] = simplified_constraints(x_next, lb, ub, data, mda_options);

    bad_constraint = max(g_next) > 1e-6 && max(g_next) > max(g_curr);
    bad_objective  = max(g_curr) <= 1e-6 && max(g_next) <= 1e-6 && f_next > f_curr;

    if bad_constraint || bad_objective
        if bad_constraint
            fprintf('Iter %d: REJECTED — constraint violated (max_g = %.4f). Shrinking move limit.\n', k, max(g_next));
        else
            fprintf('Iter %d: REJECTED — objective worsened (f = %.6f). Shrinking move limit.\n', k, f_next);
        end
        move_limit = move_limit * 0.5;
    else
        x_curr = x_next;
        fprintf('Iter %d: f = %9.6f, beta = %.4e, max_g = %9.4f, move_limit = %.3f\n', ...
            k, f_next, beta, max(g_next), move_limit);

        if k > 1 && abs(f_next - history.fval(k-1)) < conv_tol && norm(dx) < conv_tol && max(g_next) <= 1e-4
            fprintf('\nConverged.\n'); break;
        end

        if norm(dx) >= 0.99 * move_limit
            move_limit = min(0.2, move_limit * 1.2);
        end
    end
end

% =========================================================================
% POST-OPTIMIZATION: KKT & RESULTS
% =========================================================================
x_opt = denormalize_vars(x_curr, lb, ub);
f_opt = simplified_optim(x_curr, lb, ub, data, mda_options);
[g_opt, ~] = simplified_constraints(x_curr, lb, ub, data, mda_options);

fprintf('\n------------------------------------------------\n');
fprintf('OPTIMAL DESIGN VARIABLES:\n');
for i = 1:length(x_opt)
    fprintf('%-10s: %8.3f\n', var_names{i}, x_opt(i));
end
fprintf('Optimal Objective: %.6f\n', f_opt);
fprintf('------------------------------------------------\n');

[grad_f, grad_g] = get_gradients(x_curr, lb, ub, data, mda_options, h_step);

fprintf('\n--- KKT Condition Verification ---\n');
fprintf('Complementary Slackness:\n');
for j = 1:length(g_opt)
    fprintf('  %s: g = %8.4f, lambda = %8.4f\n', constraint_names{j}, g_opt(j), multipliers(j));
end

combined_grad = grad_f(:);
for j = 1:length(multipliers)
    combined_grad = combined_grad + multipliers(j) * grad_g(j, :)';
end

stationarity_error = norm(combined_grad);
fprintf('\nStationarity Error ||grad L||: %.6f\n', stationarity_error);

if stationarity_error < 0.10
    fprintf('KKT Status: Satisfied.\n');
else
    fprintf('KKT Status: Not fully satisfied.\n');
end

% =========================================================================
% PLOTS
% =========================================================================
figure('Color', 'w', 'Name', 'Optimization Convergence');
plot(history.iter, history.fval, '-bo', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on; xlabel('Iteration'); ylabel('Objective f(x)'); title('Objective Convergence');

figure('Color', 'w', 'Name', 'Design Variable Evolution', 'Position', [100, 100, 600, 800]);
for v = 1:length(x0)
    subplot(length(x0), 1, v); hold on;
    yline(lb(v), '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.2);
    yline(ub(v), '--', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.2);
    yline(x0(v), ':',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    plot(history.iter, history.x(:, v), '-o', 'LineWidth', 1.8, 'MarkerSize', 5, ...
        'Color', [0.1 0.45 0.85], 'MarkerFaceColor', [0.1 0.45 0.85]);
    yline(x_opt(v), '-', 'Color', [0.1 0.7 0.3], 'LineWidth', 1.5);
    y_range = max(history.x(:, v)) - min(history.x(:, v));
    margin  = y_range * 0.1 + 1e-6;
    ylim([min(history.x(:, v)) - margin, max(history.x(:, v)) + margin]);
    ylabel(var_names{v}); grid on;
    if v == length(x0), xlabel('Iteration'); else, set(gca, 'XTickLabel', []); end
end
sgtitle('Design Variable Evolution', 'FontWeight', 'bold', 'FontSize', 13);

figure('Color', 'w', 'Name', 'Constraint History');
if ~isempty(history.constr) && ~isempty(history.iter)
    minLen = min(length(history.iter), size(history.constr, 1));
    h = plot(history.iter(1:minLen), history.constr(1:minLen, :), '-s', 'LineWidth', 1.2);
    hold on; yline(0, 'r--', 'LineWidth', 2);
    ylim([-0.2, 0.1]); grid on;
    xlabel('Iteration'); ylabel('g(x)  [g \leq 0]'); title('Constraint Satisfaction');
    legend(h, constraint_names(1:size(history.constr, 2)), 'Location', 'bestoutside');
end

% PLOT OBJECTIVE TOPOLOGY WITH OPTIMIZATION PATH
n_points = 20; % Resolution 
v_range = linspace(lb(1), ub(1), n_points);
bpr_range = linspace(lb(2), ub(2), n_points);
[V_grid, BPR_grid] = meshgrid(v_range, bpr_range);
Z_objective = zeros(size(V_grid));
G1_clearance = zeros(size(V_grid));
G2_tit = zeros(size(V_grid));
G3_mach = zeros(size(V_grid));

for i = 1:numel(V_grid)
    x_grid = normalize_vars([V_grid(i), BPR_grid(i)], lb, ub);
    Z_objective(i) = simplified_optim(x_grid, lb, ub, data, mda_options);

    g_vals = simplified_constraints(x_grid, lb, ub, data, mda_options);
    G1_clearance(i) = g_vals(1);
    G2_tit(i) = g_vals(2);
    G3_mach(i) = g_vals(3);
end

figure('Color', 'w', 'Name', 'Design Space Topology');
hold on; grid on;
[C, h_cont] = contour(V_grid, BPR_grid, Z_objective, 15, 'k', 'LineWidth', 1.2);
clabel(C, h_cont, 'FontSize', 9, 'Color', 'k', 'LabelSpacing', 200);

[C1, h1] = contour(V_grid, BPR_grid, G1_clearance, [0 0], 'g', 'LineWidth', 2.5);
[C2, h2] = contour(V_grid, BPR_grid, G2_tit,       [0 0], 'b', 'LineWidth', 2.5);
[C3, h3] = contour(V_grid, BPR_grid, G3_mach,      [0 0], 'c', 'LineWidth', 2.5);

plot(history.x(:,1), history.x(:,2), 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', "MarkerSize", 4); % Optimization path
plot(history.x(1,1), history.x(1,2), 'gs', 'MarkerSize', 8, 'LineWidth', 2); % Start
plot(x_opt(1), x_opt(2), 'm*', 'MarkerSize', 8, 'LineWidth', 2); % Optimum

xlabel('Velocity (V)');
ylabel('Bypass Ratio (BPR)');

legend([h1, h2, h3], {'Clearance Limit (g=0)', 'TIT Limit (g=0)', 'Mach Limit (g=0)'}, ...
       'Location', 'northeastoutside');
axis([lb(1) ub(1) lb(2) ub(2)]);



















% =========================================================================
% HELPER FUNCTIONS
% =========================================================================
function x_norm = normalize_vars(x, lb, ub)
    x_norm = (x - lb) ./ (ub - lb);
end

function x_real = denormalize_vars(x_norm, lb, ub)
    x_real = x_norm .* (ub - lb) + lb;
end

function [grad_f, grad_g] = get_gradients(x, lb, ub, data, opts, h)
    n      = length(x);
    grad_f = zeros(n, 1);
    [g0, ~] = simplified_constraints(x, lb, ub, data, opts);
    grad_g  = zeros(length(g0), n);

    for i = 1:n
        x_plus  = x; x_plus(i)  = min(1, x(i) + h);
        x_minus = x; x_minus(i) = max(0, x(i) - h);

        f_plus  = simplified_optim(x_plus,  lb, ub, data, opts);
        f_minus = simplified_optim(x_minus, lb, ub, data, opts);

        actual_2h  = x_plus(i) - x_minus(i);
        grad_f(i)  = (f_plus  - f_minus)  / actual_2h;

        [g_plus,  ~] = simplified_constraints(x_plus,  lb, ub, data, opts);
        [g_minus, ~] = simplified_constraints(x_minus, lb, ub, data, opts);
        grad_g(:, i) = (g_plus - g_minus) / actual_2h;
    end
end
