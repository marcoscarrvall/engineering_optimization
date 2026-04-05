%% --- INITIAL PROBLEM INVESTIGATION ---
fprintf('Running Initial Problem Investigation...\n');
x0 = [235, 5 ]; % [V, BPR]
coefficient = 4;
lb = [x0(1) * (1 - coefficient), x0(2) * (1 - coefficient)];
ub = [x0(1) * (1 + coefficient), x0(2) * (1 + coefficient)];
options_mda.tol = 1e-8;
options_mda.max_iter = 100;
options_mda.verbose = false;
x_consts.PR_fan = 1.7;
x_consts.PR_LPC = 2.6;
x_consts.PR_HPC = 6.1;

% Define range for plotting
V_range = linspace(0, 400, 30);
BPR_range = linspace(0, 16, 30);
[VV, BB] = meshgrid(V_range, BPR_range);
F_plot = zeros(size(VV));
G1_plot = zeros(size(VV));
G2_plot = zeros(size(VV));
G3_plot = zeros(size(VV));

% Populate data for plots
for i = 1:numel(VV)
    xi = [VV(i), BB(i)];
    F_plot(i) = optim(xi, x_consts, TestAC_data, options_mda);
    [gi, ~] = constraints(xi, x_consts, TestAC_data, options_mda);
    G1_plot(i) = gi(1);
    G2_plot(i) = gi(2);
    G3_plot(i) = gi(3);
end

%% 1. Boundedness & Monotonicity (2D Contour Plot)
figure('Name', 'Design Space Investigation');
contourf(VV, BB, F_plot, 20); hold on;
colorbar; colormap jet;
title('Objective Function Contours and Constraints');
xlabel('Velocity (V)'); ylabel('Bypass Ratio (BPR)');
contour(VV, BB, G1_plot, [0 0], 'r', 'LineWidth', 2);
contour(VV, BB, G2_plot, [0 0], 'm', 'LineWidth', 2);
contour(VV, BB, G3_plot, [0 0], 'k', 'LineWidth', 2);
legend('Objective', 'g1: Clearance', 'g2: Noise', 'g3: TIT');
hold off;

%% 2. Numerical Noise & Sensitivity
V_noise_range = linspace(x0(1)-0.01, x0(1)+0.01, 100);
F_noise = zeros(1, length(V_noise_range));
for i = 1:length(V_noise_range)
    F_noise(i) = optim([V_noise_range(i), x0(2)], x_consts, TestAC_data, options_mda);
end
figure('Name', 'Numerical Noise Check');
plot(V_noise_range, F_noise, '-x');
title('Zoomed Objective Function (Noise Check)');
xlabel('Velocity (V)'); ylabel('Objective Value');
grid on;

%% 3. Sensitivity Analysis (Step Size Study)
h_sizes = logspace(-1, -10, 10);
df_dx = zeros(1, length(h_sizes));
for i = 1:length(h_sizes)
    h = h_sizes(i);
    df = (optim([x0(1)+h, x0(2)], x_consts, TestAC_data, options_mda) - optim(x0, x_consts, TestAC_data, options_mda)) / h;
    df_dx(i) = df;
end
figure('Name', 'Sensitivity Study');
loglog(h_sizes, abs(df_dx), '-o');
title('Sensitivity Study: Effect of Step Size h');
xlabel('Step Size (h)'); ylabel('|df/dV|');
grid on;

%% 4. Multi-Start Optimization with Path Plotting
starts = { [235, 11], [215, 4.5], [255, 5.5], [220, 6] };
results = zeros(length(starts), 2);
path_colors = {'w', 'c', 'y', 'g'};

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter-detailed', ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-8,...
    'MaxIterations', 500, ...          % more budget
    'FiniteDifferenceStepSize', 1e-5); % consistent FD step

% Redraw contour figure to overlay optimizer paths
figure('Name', 'Design Space Investigation - Optimizer Paths');
contourf(VV, BB, F_plot, 20); hold on;
colorbar; colormap jet;
title('Objective Function Contours, Constraints, and Optimizer Paths');
xlabel('Velocity (V)'); ylabel('Bypass Ratio (BPR)');
contour(VV, BB, G1_plot, [0 0], 'r', 'LineWidth', 2);
contour(VV, BB, G2_plot, [0 0], 'm', 'LineWidth', 2);
contour(VV, BB, G3_plot, [0 0], 'k', 'LineWidth', 2);

global path_x_global

for i = 1:length(starts)
    path_x_global = starts{i};  % reset and seed with starting point

    opts_with_path = optimoptions(options, 'OutputFcn', @trackPath);

    [x_opt, ~] = fmincon(@(x) optim(x, x_consts, TestAC_data, options_mda), starts{i}, ...
                         [], [], [], [], lb, ub, ...
                         @(x) constraints(x, x_consts, TestAC_data, options_mda), ...
                         opts_with_path);
    results(i, :) = x_opt;

    path_x = path_x_global;

    % Plot path
    col = path_colors{i};
    plot(path_x(:,1), path_x(:,2), '-o', 'Color', col, 'LineWidth', 1.5, ...
         'MarkerSize', 4, 'MarkerFaceColor', col);
    % Mark start (square) and optimum (pentagram)
    plot(starts{i}(1), starts{i}(2), 's', 'Color', col, 'MarkerSize', 9, ...
         'MarkerFaceColor', col, 'LineWidth', 1.5);
    plot(x_opt(1), x_opt(2), 'p', 'Color', col, 'MarkerSize', 12, ...
         'MarkerFaceColor', col, 'LineWidth', 1.5);
end

legend('Objective', 'g1: Clearance', 'g2: Noise', 'g3: TIT', ...
       'Path 1', 'Path 2', 'Path 3', 'Path 4', 'Location', 'best');
hold off;

% Check if the results are the same
disp('Comparison of optima from different starts:');
disp(results);

function stop = trackPath(x, ~, ~)
    global path_x_global
    path_x_global = [path_x_global; x];
    stop = false;
end