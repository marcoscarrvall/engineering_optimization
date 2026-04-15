clear; clc; close all;

% --- Design Variables: [V, BPR, PR_fan, PR_LPC, PR_HPC] ---
x0 = [230, 8.0, 1.55, 1.55, 22.0];
lb = [100, 1,   0.5,  0.5,  4  ];
ub = [300, 30,  3.0,  4.0,  35 ];

x0_norm = normalize_vars(x0, lb, ub);
lb_norm = zeros(size(x0_norm));
ub_norm = ones(size(x0_norm));

mda_options.tol      = 1e-3;
mda_options.max_iter = 100;
mda_options.verbose  = false;

data = A320data();

var_names = {'V', 'BPR', 'PR_{fan}', 'PR_{LPC}', 'PR_{HPC}'};
units     = {'[m/s]', '[-]', '[-]', '[-]', '[-]'};
con_names = {'Ground Clearance', 'Max TIT', 'Fan Tip Mach'};

% =========================================================================
% 1. BOUNDEDNESS SWEEP (one variable at a time)
% =========================================================================
f0       = optim(x0_norm, lb, ub, data, mda_options);
[g0, ~]  = constraints(x0_norm, lb, ub, data, mda_options);

figure('Name', 'Boundedness Sweep');
for i = 1:5
    x_test_norm = linspace(0, 1, 50);
    x_test_phys = x_test_norm .* (ub(i) - lb(i)) + lb(i);

    f_vals = zeros(1, 50);
    g_vals = zeros(50, 3);

    for j = 1:50
        x_eval    = x0_norm;
        x_eval(i) = x_test_norm(j);
        f_vals(j) = optim(x_eval, lb, ub, data, mda_options);
        [g, ~]    = constraints(x_eval, lb, ub, data, mda_options);
        g_vals(j, :) = g';
    end

    subplot(2, 3, i);

    yyaxis left
    plot(x_test_phys, f_vals, 'b-', 'LineWidth', 1.5); hold on;
    plot(x0(i), f0, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    ylabel('Objective f');
    xlim([lb(i), ub(i)]);

    yyaxis right
    p1 = plot(x_test_phys, g_vals(:,1), 'r--', 'LineWidth', 1.0); hold on;
    p2 = plot(x_test_phys, g_vals(:,2), 'g--', 'LineWidth', 1.0);
    p3 = plot(x_test_phys, g_vals(:,3), 'm--', 'LineWidth', 1.0);
    plot(x0(i), g0(1), 'ks', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
    plot(x0(i), g0(2), 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'HandleVisibility', 'off');
    plot(x0(i), g0(3), 'ks', 'MarkerFaceColor', 'm', 'MarkerSize', 6, 'HandleVisibility', 'off');
    yline(0, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    ylim([-0.3, 0.2]);
    ylabel('Constraint g');

    title(['Boundedness: ', var_names{i}]);
    xlabel([var_names{i}, ' ', units{i}]);
    grid on;
    legend([p1, p2, p3], con_names, 'Location', 'northeast');
end

% =========================================================================
% 2. CONTOUR PLOT: Objective & Feasible Boundary (V vs BPR)
% =========================================================================
V_range   = linspace(lb(1), ub(1), 30);
BPR_range = linspace(lb(2), ub(2), 30);
[V_grid, BPR_grid] = meshgrid(V_range, BPR_range);

F_grid     = zeros(size(V_grid));
G_max_grid = zeros(size(V_grid));

for r = 1:size(V_grid, 1)
    for c = 1:size(V_grid, 2)
        x_phys    = x0;
        x_phys(1) = V_grid(r, c);
        x_phys(2) = BPR_grid(r, c);
        x_norm_rc = normalize_vars(x_phys, lb, ub);

        F_grid(r, c)     = optim(x_norm_rc, lb, ub, data, mda_options);
        [g, ~]           = constraints(x_norm_rc, lb, ub, data, mda_options);
        G_max_grid(r, c) = max(g);
    end
end

figure('Name', 'Objective Contours with Feasible Boundary');
contour(V_grid, BPR_grid, G_max_grid, [0 0], 'r', 'LineWidth', 2); hold on;
contour(V_grid, BPR_grid, F_grid, 20, 'LineWidth', 1.2);
cb = colorbar('westoutside');
cb.Label.String = 'Objective Value f';
plot(x0(1), x0(2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
xlabel('Velocity V [m/s]');
ylabel('Bypass Ratio BPR [-]');
title('Objective Contours with Feasible Boundary');
legend('Feasible Boundary (g = 0)', 'Objective Contours', 'Initial Point', 'Location', 'northeast');
grid on;

% =========================================================================
% 3. NUMERICAL NOISE SCAN
% =========================================================================
lb = [200, 5.0, 1.1, 1.1, 8.0];
ub = [235, 14.0, 1.7, 2.0, 25.0];

n_points = 1000;

figure('Name', 'Numerical Noise Scan');
for i = 1:5
    x_test_norm = linspace(0, 1, n_points);
    f_noise     = zeros(1, n_points);

    for j = 1:n_points
        x_eval    = x0_norm;
        x_eval(i) = x_test_norm(j);
        f_noise(j) = optim(x_eval, lb, ub, data, mda_options);
    end

    x_phys_plot = x_test_norm .* (ub(i) - lb(i)) + lb(i);

    subplot(2, 3, i);
    plot(x_phys_plot, f_noise, 'b-');
    title(['Noise Scan: ', var_names{i}]);
    xlabel([var_names{i}, ' ', units{i}]);
    xlim([lb(i), ub(i)]);
    ylabel('Objective f');
    grid on;
end

% =========================================================================
% 4. FINITE DIFFERENCE STEP SIZE STUDY
% =========================================================================
h_range  = logspace(-1, -11, 40);
n_vars   = length(var_names);
n_h      = length(h_range);

dfdx_central = zeros(n_vars, n_h);
dfdx_forward = zeros(n_vars, n_h);

figure('Name', 'Step Size Sensitivity Study', 'Position', [100, 100, 1200, 800]);
for i = 1:n_vars
    f0 = optim(x0_norm, lb, ub, data, mda_options);

    for j = 1:n_h
        h = h_range(j);

        x_plus    = x0_norm; x_plus(i)  = min(1, x0_norm(i) + h);
        x_minus   = x0_norm; x_minus(i) = max(0, x0_norm(i) - h);
        f_plus    = optim(x_plus,  lb, ub, data, mda_options);
        f_minus   = optim(x_minus, lb, ub, data, mda_options);

        dfdx_forward(i, j) = (f_plus - f0)           / h;
        dfdx_central(i, j) = (f_plus - f_minus)       / (2 * h);
    end

    subplot(2, 3, i);
    loglog(h_range, abs(dfdx_forward(i, :)), 'r--', 'LineWidth', 1.1); hold on;
    loglog(h_range, abs(dfdx_central(i, :)), 'b-',  'LineWidth', 1.5);
    xline(1e-3, 'k:',  'MDA Tol',  'LabelVerticalAlignment', 'bottom');
    xline(5e-3, 'g--', 'Chosen h', 'LabelVerticalAlignment', 'bottom');
    grid on;
    title(['Step Size: ', var_names{i}]);
    xlabel('Step size h');
    ylabel(['|df/d', var_names{i}, '|']);
    legend('Forward', 'Central', 'Location', 'best');
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function x_norm = normalize_vars(x, lb, ub)
    x_norm = (x - lb) ./ (ub - lb);
end

function x_real = denormalize_vars(x_norm, lb, ub)
    x_real = x_norm .* (ub - lb) + lb;
end