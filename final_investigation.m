clear; clc; close all;

constraint_names = {'clearance', 'noise', 'TIT'};

% Initial guess and bounds
x0 = [235, 11, 1.7, 2.6, 6.1]; 

options_mda.tol = 1e-6;
options_mda.max_iter = 100;
options_mda.verbose = false;

coefficient = 1;
lb = x0 * (1 - coefficient);
ub = x0 * (1 + coefficient);

% Global storage
global optHistory current_c
optHistory.fval  = [];
optHistory.constr = [];
optHistory.iter  = [];
optHistory.x     = [];   % store full x at each iteration for path plotting

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-8,...
    'MaxFunctionEvaluations', 5000, ...
    'FiniteDifferenceStepSize', 1e-8, ...
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

% --- 2D DESIGN SPACE PLOTS (V vs each other variable) ---
% Others fixed at x0 (baseline). Optimizer path overlaid on each panel.
% --- 2D DESIGN SPACE PLOTS (V vs each other variable) ---
N = 25;

var_labels = {'BPR', 'PR_{fan}', 'PR_{LPC}', 'PR_{HPC}'};
var_idx    = [2, 3, 4, 5];

V_range = linspace(lb(1), ub(1), N);
path_V  = optHistory.x(:, 1);

figure('Color', 'w', 'Name', 'Design Space: V vs Each Variable');

for p = 1:4
    vi      = var_idx(p);
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

    % --- Sanitize: strip complex parts produced by MDA divergence ---
    F_surf  = real(F_surf);
    G1_surf = real(G1_surf);
    G2_surf = real(G2_surf);
    G3_surf = real(G3_surf);

    % Mask extreme outliers so the colorbar is not dominated by bad points
    F_med = median(F_surf(:), 'omitnan');
    F_std = std(F_surf(:),    'omitnan');
    F_surf(abs(F_surf - F_med) > 10 * F_std) = NaN;

    subplot(2, 2, p);
    hold on;

    contourf(VV, YY, F_surf, 20, 'LineStyle', 'none');
    colormap(gca, jet);
    cb = colorbar;
    cb.Label.String = 'f(x)';

    [~, hg1] = contour(VV, YY, G1_surf, [0 0], 'r', 'LineWidth', 2);
    [~, hg2] = contour(VV, YY, G2_surf, [0 0], 'm', 'LineWidth', 2);
    [~, hg3] = contour(VV, YY, G3_surf, [0 0], 'k', 'LineWidth', 2);

    path_y = optHistory.x(:, vi);
    h_path = plot(path_V, path_y, 'w-o', ...
                  'LineWidth', 1.8, 'MarkerSize', 5, 'MarkerFaceColor', 'w');

    h_start = plot(x0(1), x0(vi), 'g^', ...
                   'MarkerSize', 11, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);

    h_opt = plot(x_opt(1), x_opt(vi), 'wp', ...
                 'MarkerSize', 15, 'MarkerFaceColor', 'y', 'LineWidth', 1.5);

    xlabel('V (m/s)');
    ylabel(var_labels{p});
    title(['V vs ' var_labels{p}]);
    legend([hg1, hg2, hg3, h_path, h_start, h_opt], ...
           'g1: clearance', 'g2: noise', 'g3: TIT', ...
           'Optimizer path', 'x_0', 'x_{opt}', ...
           'Location', 'best', 'TextColor', 'w');
    grid on;
    hold off;
end

sgtitle('Design Space: V vs Each Variable (others fixed at x_0)', ...
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
        optHistory.x      = [optHistory.x;       x(:)'];  % log full x vector
    end
end