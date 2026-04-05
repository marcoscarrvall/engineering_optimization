clear; clc; close all;

constraint_names = {'clearance', 'noise', 'TIT'};

% Initial guess and bounds
x0 = [235, 9.6, 1.55, 1.55, 22.2]; 

options_mda.tol = 1e-8;
options_mda.max_iter = 200;
options_mda.verbose = false;

coefficient = 1;
lb = [50, 4, 1.1, 1.1, 10];
ub = [260, 20, 1.8, 1.8, 25];

% Global storage — all fields including x
global optHistory current_c
optHistory.fval   = [];
optHistory.constr = [];
optHistory.iter   = [];
optHistory.x      = [];

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...            % SQP is good for physical constraints
    'Display', 'iter-detailed', ...    % Gives more info on which constraint is the "killer"
    'OptimalityTolerance', 1e-6, ...   % Relaxed
    'StepTolerance', 1e-6, ...         % Relaxed
    'ConstraintTolerance', 1e-4, ...   % Allows a tiny bit of "breathing room"
    'MaxFunctionEvaluations', 10000, ...
    'OutputFcn', @outfun);

fprintf('\n%-10s %-12s %-12s %-12s\n', 'Iter', 'V (x1)', 'BPR (x2)', 'f(x)');
fprintf('-----------------------------------------------------------\n');

[x_opt, f_opt] = fmincon(@(x) optim(x, TestAC_data, options_mda), x0, ...
                         [], [], [], [], lb, ub, ...
                         @(x) constraints_catcher(x, TestAC_data, options_mda), options);

fprintf('\nOptimal: V=%.3f  BPR=%.3f  PR_fan=%.3f  PR_LPC=%.3f  PR_HPC=%.3f  f=%.6f\n', ...
    x_opt(1), x_opt(2), x_opt(3), x_opt(4), x_opt(5), f_opt);

% --- CONVERGENCE PLOTS ---
figure('Color', 'w', 'Name', 'Optimization Convergence');

subplot(2,1,1);
plot(optHistory.iter, optHistory.fval, '-bo', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on; ylabel('Objective f(x)'); title('Minimize Objective');

subplot(2,1,2);
if ~isempty(optHistory.constr)
    h = plot(optHistory.iter, optHistory.constr, '-s', 'LineWidth', 1.2);
    hold on;
    yline(0, 'r--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
    grid on; ylabel('Value (c \leq 0)'); xlabel('Iteration');
    title('Constraint Satisfaction');
    num_found = size(optHistory.constr, 2);
    if length(constraint_names) >= num_found
        legend(h, constraint_names(1:num_found), 'Location', 'bestoutside');
    else
        legend(h, 'Location', 'bestoutside');
    end
end

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

% --- 2D DESIGN SPACE PLOTS WITH FIXED VARIABLES ---
N = 25;
var_names_plot = {'BPR', 'PR_{fan}', 'PR_{LPC}', 'PR_{HPC}'};
var_idx        = [2, 3, 4, 5];

figure('Color', 'w', 'Name', 'Design Space: V vs Each Variable (2D subspace)');

for p = 1:4
    vi = var_idx(p);

    % --- Build grid for background contour ---
    V_range = linspace(lb(1), ub(1), N);
    y_range = linspace(lb(vi), ub(vi), N);
    [VV, YY] = meshgrid(V_range, y_range);

    F_surf  = zeros(N, N);
    G1_surf = zeros(N, N);
    G2_surf = zeros(N, N);
    G3_surf = zeros(N, N);

    for r = 1:N
        for c = 1:N
            xi     = x0;
            xi(1)  = VV(r, c);
            xi(vi) = YY(r, c);

            F_surf(r,c)  = optim(xi, TestAC_data, options_mda);
            [gi, ~]      = constraints(xi, TestAC_data, options_mda);
            G1_surf(r,c) = gi(1);
            G2_surf(r,c) = gi(2);
            G3_surf(r,c) = gi(3);
        end
    end

    % Sanitize complex/diverged values
    F_surf  = real(F_surf);
    G1_surf = real(G1_surf);
    G2_surf = real(G2_surf);
    G3_surf = real(G3_surf);
    F_med = median(F_surf(:), 'omitnan');
    F_std = std(F_surf(:),    'omitnan');
    F_surf(abs(F_surf - F_med) > 10 * F_std) = NaN;

    % Normalize F_surf to [-1, 1]
    F_min  = min(F_surf(:), [], 'omitnan');
    F_max  = max(F_surf(:), [], 'omitnan');
    F_norm = 2 * (F_surf - F_min) / (F_max - F_min) - 1;

    % --- Run 2D-constrained optimization ---
    lb_2d = lb;
    ub_2d = ub;
    for k = 1:length(x0)
        if k ~= 1 && k ~= vi
            lb_2d(k) = x0(k);
            ub_2d(k) = x0(k);
        end
    end

    global path2d_x
    path2d_x = x0;

    opts_2d = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'none', ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-8, ...
        'MaxFunctionEvaluations', 5000, ...
        'FiniteDifferenceStepSize', 1e-8, ...
        'OutputFcn', @outfun_2d);

    [x_opt_2d, ~] = fmincon( ...
        @(x) optim(x, TestAC_data, options_mda), x0, ...
        [], [], [], [], lb_2d, ub_2d, ...
        @(x) constraints(x, TestAC_data, options_mda), opts_2d);

    path2d = path2d_x;

    % --- Plot ---
    subplot(2, 2, p);
    hold on;

    contourf(VV, YY, F_norm, 20, 'LineStyle', 'none');
    colormap(gca, jet);
    cb = colorbar;
    cb.Label.String = 'f(x) normalised [-1, 1]';
    clim([-1, 1]);

    % Constraint boundaries
    legend_handles = [];
    legend_labels  = {};
    try
        [~, hg1] = contour(VV, YY, G1_surf, [0 0], 'r', 'LineWidth', 2);
        legend_handles(end+1) = hg1; legend_labels{end+1} = 'g1: clearance';
    catch; end
    try
        [~, hg2] = contour(VV, YY, G2_surf, [0 0], 'm', 'LineWidth', 2);
        legend_handles(end+1) = hg2; legend_labels{end+1} = 'g2: noise';
    catch; end
    try
        [~, hg3] = contour(VV, YY, G3_surf, [0 0], 'k', 'LineWidth', 2);
        legend_handles(end+1) = hg3; legend_labels{end+1} = 'g3: TIT';
    catch; end

    % Optimizer path
    h_path = plot(path2d(:,1), path2d(:,vi), 'w-', 'LineWidth', 2);
    scatter(path2d(:,1), path2d(:,vi), 40, (1:size(path2d,1))', ...
            'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);

    % Start and optimum markers
    h_start = plot(x0(1), x0(vi), 'g^', ...
                   'MarkerSize', 11, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
    h_opt   = plot(x_opt_2d(1), x_opt_2d(vi), 'yp', ...
                   'MarkerSize', 15, 'MarkerFaceColor', 'y', 'LineWidth', 1.5);

    % Axis limits locked to bounds
    xlim([lb(1),  ub(1)]);
    ylim([lb(vi), ub(vi)]);

    legend_handles = [legend_handles, h_path, h_start, h_opt];
    legend_labels  = [legend_labels, 'Optimizer path', 'x_0', 'x_{opt}'];
    legend(legend_handles, legend_labels, ...
           'Location', 'best', 'TextColor', 'w', 'Color', [0.15 0.15 0.15]);

    xlabel('V (m/s)');
    ylabel(var_names_plot{p});
    title(['V vs ' var_names_plot{p} ' (others fixed at x_0)']);
    grid on;
    hold off;
end

sgtitle('2D Design Space: Constrained Optimizer Path in Each Subspace', ...
        'FontWeight', 'bold', 'FontSize', 13);

% --- INTERNAL FUNCTIONS ---

function [c, ceq] = constraints_catcher(x, data, opts)
    global current_c
    [c, ceq] = constraints(x, data, opts);
    current_c = c;
end

function stop = outfun(x, optimValues, state)
    stop = false;
    global optHistory current_c
    if strcmp(state, 'iter')
        fprintf('%-10d %-12.4f %-12.4f %-12.6f\n', ...
            optimValues.iteration, x(1), x(2), optimValues.fval);
        optHistory.iter   = [optHistory.iter;   optimValues.iteration];
        optHistory.fval   = [optHistory.fval;   optimValues.fval];
        optHistory.constr = [optHistory.constr;  current_c(:)'];
        optHistory.x      = [optHistory.x;       x(:)'];
    end
end

function stop = outfun_2d(x, optimValues, state)
    stop = false;
    global path2d_x
    if strcmp(state, 'iter')
        path2d_x = [path2d_x; x(:)'];
    end
end
