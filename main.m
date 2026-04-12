clear; clc; close all;

% Initial Vector [V, BPR, PR_fan, PR_LPC, PR_HPC]
x0 = [230, 8.0, 1.55, 1.55, 22.0];

lb = [200, 5.0, 1.1, 1.1, 8];
ub = [235, 14.0, 1.7, 2, 25.0];


x0_norm = normalize_vars(x0, lb, ub);

lb_norm = zeros(size(x0_norm));
ub_norm = ones(size(x0_norm));


% No constraints for now
A = []; b = []; Aeq = []; beq = [];

optimizer_options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ... % More stable than SQP for complex MDA
    'Display', 'iter-detailed', ...    % Gives you more info on why steps are taken
    'FiniteDifferenceType', 'central',... % Higher accuracy for gradients
    'FiniteDifferenceStepSize', 1e-2, ... % Smaller, but still above MDA noise
    'OptimalityTolerance', 1e-6, ...    % Don't chase "ghost" precision
    'ConstraintTolerance', 1e-5, ...    % Level of acceptable constraint violation
    'StepTolerance', 1e-6, ...          % Stop if x barely changes
    'UseParallel', false);              % KEEP FALSE if using global optHistory

global optHistory
optHistory.iter   = [];
optHistory.fval   = [];
optHistory.x      = [];
optHistory.constr = [];

optHistory.count  = 0;

% MDA options
mda_options.tol = 1e-3;
mda_options.max_iter = 100;
mda_options.verbose = false;

% Load data
data = A320data;

% Run optimization
[x_opt, f_opt] = fmincon(@(x) optim(x, lb, ub, data, mda_options), x0_norm, A, b, Aeq, beq, lb_norm, ub_norm,@(x) constraints(x, lb, ub, data, mda_options), optimizer_options);

function x_real = denormalize_vars(x_norm, lb, ub)
    % Scales [0, 1] values back to physical units (for the MDA)
    x_real = x_norm .* (ub - lb) + lb;
end
x_opt = denormalize_vars(x_opt, lb, ub);

fprintf('\nOptimal Design Variables:\n');
var_names = {'V (m/s)', 'BPR', 'PR_{fan}', 'PR_{LPC}', 'PR_{HPC}'};
for i = 1:length(x_opt)
    fprintf('%s: %.3f\n', var_names{i}, x_opt(i));
end

fprintf('Optimal Objective (Relative Range Error): %.6f\n', f_opt);
% --- CONVERGENCE PLOTS ---
figure('Color', 'w', 'Name', 'Optimization Convergence');

hold on; % Allow multiple items on one plot

% Plot the history
plot(optHistory.iter, optHistory.fval, '-bo', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlim([0, max(optHistory.iter)]);
grid on; 
ylabel('Objective f(x)'); 
title('Minimize Objective');


% --- DESIGN VARIABLE EVOLUTION ---
var_names = {'V (m/s)', 'BPR', 'PR_{fan}', 'PR_{LPC}', 'PR_{HPC}'};
n_vars    = length(x0);

figure('Color', 'w', 'Name', 'Design Variable Evolution');

for v = 1:n_vars
    subplot(n_vars, 1, v);
    hold on;

    yline(lb(v), '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.2, ...
          'Label', 'lb', 'LabelVerticalAlignment', 'bottom');
    yline(ub(v), '--', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.2, ...
          'Label', 'ub', 'LabelVerticalAlignment', 'top');
    yline(x0(v), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, ...
          'Label', 'x_0', 'LabelVerticalAlignment', 'middle');

    plot(optHistory.iter, optHistory.x(:, v), '-o', ...
         'LineWidth', 1.8, 'MarkerSize', 5, ...
         'Color', [0.1 0.45 0.85], 'MarkerFaceColor', [0.1 0.45 0.85]);

    yline(x_opt(v), '-', 'Color', [0.1 0.7 0.3], 'LineWidth', 1.5, ...
          'Label', 'x_{opt}', 'LabelVerticalAlignment', 'middle');

    y_min = min(optHistory.x(:, v));
    y_max = max(optHistory.x(:, v));
    margin = (y_max - y_min) * 0.1;
    if margin == 0
        margin = abs(y_min) * 0.05 + 1e-6;
    end
    ylim([y_min - margin, y_max + margin]);

    ylabel(var_names{v});
    grid on;
    if v == n_vars
        xlabel('Iteration');
    else
        set(gca, 'XTickLabel', []);
    end
end

sgtitle('Design Variable Evolution per Iteration', ...
        'FontWeight', 'bold', 'FontSize', 13);


constraint_names = {'Clearance', 'TIT Limit', 'Tip Mach'};

figure('Color', 'w', 'Name', 'Optimization Constraints');

if ~isempty(optHistory.constr) && ~isempty(optHistory.iter)
    % Find the minimum length between the two to prevent the error
    minLen = min(length(optHistory.iter), size(optHistory.constr, 1));
    
    % Synchronize data
    x_plot = optHistory.iter(1:minLen);
    y_plot = optHistory.constr(1:minLen, :);
    
    h = plot(x_plot, y_plot, '-s', 'LineWidth', 1.2);
    ylim([-1, 1]);
    hold on;
    yline(0, 'r--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
    grid on; 
    ylabel('Value (c \leq 0)'); 
    xlabel('Evaluation Count');
    title('Constraint Satisfaction');
    
    % Legend logic
    num_found = size(y_plot, 2);
    if exist('constraint_names', 'var') && length(constraint_names) >= num_found
        legend(h, constraint_names(1:num_found), 'Location', 'bestoutside');
    end
end



function x_norm = normalize_vars(x, lb, ub)
    % Scales real values to the [0, 1] range
    x_norm = (x - lb) ./ (ub - lb);
end


% --- FIND DESIGNS WHERE FVAL IS 1 ---
% Using a small tolerance (1e-4) in case of numerical precision
target_val = 1;
tol = 1e-4;

% Find indices where fval is approximately 1
idx_fail = find(abs(optHistory.fval - target_val) < tol);

if ~isempty(idx_fail)
    fprintf('\nFound %d design vectors where fval is approx 1 (Potential MDA Failures):\n', length(idx_fail));
    fprintf('----------------------------------------------------------------------\n');
    fprintf('Iter |   V   |  BPR  | PR_fan | PR_LPC | PR_HPC | Objective\n');
    fprintf('----------------------------------------------------------------------\n');
    for i = 1:length(idx_fail)
        idx = idx_fail(i);
        x_fail = optHistory.x(idx, :);
        % Display physical values
        fprintf('%4d | %5.1f | %5.1f | %6.3f | %6.3f | %6.1f | %9.4f\n', ...
            idx, x_fail(1), x_fail(2), x_fail(3), x_fail(4), x_fail(5), optHistory.fval(idx));
    end
else
    fprintf('\nNo design vectors found with fval = 1.\n');
end