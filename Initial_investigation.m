%% --- INITIAL PROBLEM INVESTIGATION ---
fprintf('Running Initial Problem Investigation...\n');

x0 = [235, 5 ]; % [V, BPR]

coefficient = 4;

lb = [x0(1) * (1 - coefficient), x0(2) * (1 - coefficient)];
ub = [x0(1) * (1 + coefficient), x0(2) * (1 + coefficient)];

options_mda.tol = 1e-6;
options_mda.max_iter = 100;
options_mda.verbose = false;

x_consts.PR_fan = 1.7;
x_consts.PR_LPC = 2.6;
x_consts.PR_HPC = 6.1;

% Define range for plotting (10% around x0 as per your code)
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

% Overlay constraints (Negative Null Form: g <= 0 is feasible)
% We plot the 0-level line where the constraint becomes active
contour(VV, BB, G1_plot, [0 0], 'r', 'LineWidth', 2); % Constraint 1
contour(VV, BB, G2_plot, [0 0], 'm', 'LineWidth', 2); % Constraint 2
contour(VV, BB, G3_plot, [0 0], 'k', 'LineWidth', 2); % Constraint 3
legend('Objective', 'g1: Clearance', 'g2: Noise', 'g3: TIT');
hold off;

%% 2. Numerical Noise & Sensitivity
% Pick one variable (Velocity) and zoom in significantly to see MDA noise
% Preallocate with zeros based on the size of your range
V_noise_range = linspace(x0(1)-0.01, x0(1)+0.01, 100);
F_noise = zeros(1, length(V_noise_range)); 

for i = 1:length(V_noise_range)
    % Use the index 'i' to fill the preallocated spots
    F_noise(i) = optim([V_noise_range(i), x0(2)], x_consts, TestAC_data, options_mda);
end

figure('Name', 'Numerical Noise Check');
plot(V_noise_range, F_noise, '-x');
title('Zoomed Objective Function (Noise Check)');
xlabel('Velocity (V)'); ylabel('Objective Value');
grid on;
% Note: If this line is jagged, your MDA tolerance is too loose.

%% 3. Sensitivity Analysis (Step Size Study)
h_sizes = logspace(-1, -10, 10);
df_dx = zeros(1, length(h_sizes)); % Preallocate

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

% Define different starting points
starts = { [235, 5], [215, 4.5], [255, 5.5], [220, 6] };
results = zeros(length(starts), 2);

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...        
    'Display', 'iter-detailed', ... 
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6);

    
for i = 1:length(starts)
    [x_opt, f_opt] = fmincon(@(x) optim(x, x_consts, TestAC_data, options_mda), starts{i}, ...
                             [], [], [], [], lb, ub, @(x) constraints(x, x_consts, TestAC_data, options_mda), options);
    results(i, :) = x_opt;
end

% Check if the results are the same
disp('Comparison of optima from different starts:');
disp(results);
