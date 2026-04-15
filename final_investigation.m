clear; clc; close all;

% Initial Vector [V, BPR, PR_fan, PR_LPC, PR_HPC]
x0 = [230, 8.0, 1.55, 1.55, 22.0];

lb = [100, 1, 0.5, 0.5, 4];
ub = [300, 30, 3, 4, 35];

x0_norm = normalize_vars(x0, lb, ub);

lb_norm = zeros(size(x0_norm));
ub_norm = ones(size(x0_norm));

% No constraints for now
A = []; b = []; Aeq = []; beq = [];

optimizer_options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ... 
    'Display', 'iter-detailed', ...    
    'FiniteDifferenceType', 'central',... 
    'FiniteDifferenceStepSize', 1e-2, ... 
    'OptimalityTolerance', 1e-6, ...    
    'ConstraintTolerance', 1e-5, ...    
    'StepTolerance', 1e-6, ...          
    'UseParallel', false);              

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

% Initial setup
var_names = {'V', 'BPR', 'PR_{fan}', 'PR_{LPC}', 'PR_{HPC}'};
units = {'[m/s]', '[-]', '[-]', '[-]', '[-]'};
con_names = {'Ground Clearance', 'Max TIT', 'Fan Tip Mach'};

% Pre-calculate initial objective and constraints for the reference point
f0 = optim(x0_norm, lb, ub, data, mda_options);
[g0, ~] = constraints(x0_norm, lb, ub, data, mda_options);

figure;
for i = 1:5
    x_test_norm = linspace(0, 1, 50); 
    
    % Denormalize for plotting purposes (X-axis)
    x_test_phys = x_test_norm .* (ub(i) - lb(i)) + lb(i);
    
    f_vals = zeros(1, 50);
    g_vals = zeros(50, 3); 
    
    for j = 1:50
        x_eval = x0_norm; 
        x_eval(i) = x_test_norm(j); 
        f_vals(j) = optim(x_eval, lb, ub, data, mda_options); 
        [g, ~] = constraints(x_eval, lb, ub, data, mda_options);
        g_vals(j, :) = g'; 
    end
    
    subplot(2, 3, i);
    
    % --- Plot Objective (Left Axis) ---
    yyaxis left
    plot(x_test_phys, f_vals, 'b-', 'LineWidth', 1.5); hold on;
    % Plot the Initial Point on the Objective curve at its physical location
    plot(x0(i), f0, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Initial Point (f)');
    ylabel('Objective (f)');
    xlim([lb(i) ub(i)]);
    % --- Plot Constraints (Right Axis) ---
    yyaxis right
    p1 = plot(x_test_phys, g_vals(:,1), 'r--', 'LineWidth', 1.0); hold on;
    p2 = plot(x_test_phys, g_vals(:,2), 'g--', 'LineWidth', 1.0);
    p3 = plot(x_test_phys, g_vals(:,3), 'm--', 'LineWidth', 1.0);
    
    % Plot the Initial Point for each constraint at its physical location
    plot(x0(i), g0(1), 'ks', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
    plot(x0(i), g0(2), 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'HandleVisibility', 'off');
    plot(x0(i), g0(3), 'ks', 'MarkerFaceColor', 'm', 'MarkerSize', 6, 'HandleVisibility', 'off');

    % Add horizontal boundary line (Negative Null Form: g <= 0)
    line([lb(i) ub(i)], [0 0], 'Color', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    ylim([-0.3, 0.2]); 
    ylabel('Constraint Value (g)');
    
    title(['Boundedness: ', var_names{i}]);
    xlabel([var_names{i}, ' ', units{i}]);
    grid on;
    legend([p1, p2, p3], con_names, 'Location', 'northeast');
end

function x_norm = normalize_vars(x, lb, ub)
    % Scales real values to the [0, 1] range
    x_norm = (x - lb) ./ (ub - lb);
end
function x_real = denormalize_vars(x_norm, lb, ub)
    % Scales [0, 1] values back to physical units (for the MDA)
    x_real = x_norm .* (ub - lb) + lb;
end
% Define grid for V and BPR
V_range = linspace(lb(1), ub(1), 30);
BPR_range = linspace(lb(2), ub(2), 30);
[V_grid, BPR_grid] = meshgrid(V_range, BPR_range);

F_grid = zeros(size(V_grid));
G_max_grid = zeros(size(V_grid));

for r = 1:size(V_grid, 1)
    for c = 1:size(V_grid, 2)
        % Set current V and BPR, others remain at x0_norm baseline
        x_phys = x0; 
        x_phys(1) = V_grid(r,c);
        x_phys(2) = BPR_grid(r,c);
        
        x_norm = (x_phys - lb) ./ (ub - lb);
        
        F_grid(r,c) = optim(x_norm, lb, ub, data, mda_options);
        [g, ~] = constraints(x_norm, lb, ub, data, mda_options);
        G_max_grid(r,c) = max(g); % Highest violation
    end
end

% Plotting
figure;
% Overlay the feasible boundary (where max constraint is 0)
contour(V_grid, BPR_grid, G_max_grid, [0 0], 'r', 'LineWidth', 2); 

hold on;
% --- Plot Objective (Left Axis) ---
yyaxis left
% 'ShowText', 'off' removes the numbers from the lines
[C, h_obj_cont] = contour(V_grid, BPR_grid, F_grid, 20, 'ShowText', 'off', 'LineWidth', 1.2); 
hold on;

% Add a colorbar to show the scale of the objective function
cb = colorbar('westoutside'); % Places the scale on the far left
cb.Label.String = 'Objective Value (f)';

% Plot the Initial Point
plot(x0(1), x0(2), 'ko', 'MarkerFaceColor', 'white', 'MarkerSize', 8);
ylabel('Velocity V [m/s]');
xlabel('Bypass Ratio (BPR) [-]');

legend('Feasible Boundary (g=0)', 'Objective Contours', 'Initial Point', 'Location', 'northeast');
title('Objective Function Contours with Feasible Boundary');



% Setup for the scan
var_names = {'V', 'BPR', 'PR_{fan}', 'PR_{LPC}', 'PR_{HPC}'};
delta = 1; % Tiny range to look for noise
n_points = 1000;
noise_fig = figure('Name', 'Numerical Noise Investigation');
lb = [200, 5.0, 1.1, 1.1, 8];
ub = [235, 14.0, 1.7, 2, 25.0];

for i = 1:5
    % Define the tiny interval around the normalized initial point
    x_test_norm = linspace(max(0, x0_norm(i) - delta), min(1, x0_norm(i) + delta), n_points);
    f_noise = zeros(1, n_points);
    
    for j = 1:n_points
        x_eval = x0_norm; 
        x_eval(i) = x_test_norm(j); 
        f_noise(j) = optim(x_eval, lb, ub, data, mda_options); 
    end
    
    subplot(2, 3, i);
    
    % FIX: Only denormalize the i-th variable being plotted
    x_phys_plot = x_test_norm .* (ub(i) - lb(i)) + lb(i);
    
    plot(x_phys_plot, f_noise, 'b-'); % Removed '.' to see the line smoothness better
    title(['Noise Scan: ', var_names{i}]);
    xlabel(['Physical Value ', units{i}]); % Use physical units for the report
    xlim([lb(i) ub(i)]);
    ylabel('Objective f');
    grid on;
end

% --- Setup and Pre-allocation ---
var_names = {'V', 'BPR', 'PR_{fan}', 'PR_{LPC}', 'PR_{HPC}'};
h_range = logspace(-1, -11, 40); 
n_vars = length(var_names);
n_h = length(h_range);

% Matrix to store results for ALL variables
all_dfdx_central = zeros(n_vars, n_h); 
all_dfdx_forward = zeros(n_vars, n_h);

figure('Name', 'Sensitivity Step Size Study', 'Position', [100, 100, 1200, 800]);

% --- FIRST LOOP: Data Collection ---
for i = 1:n_vars
    f0 = optim(x0_norm, lb, ub, data, mda_options);
    
    for j = 1:n_h
        h = h_range(j);
        
        % Forward
        x_plus = x0_norm; 
        x_plus(i) = min(1, x0_norm(i) + h); 
        f_plus = optim(x_plus, lb, ub, data, mda_options);
        all_dfdx_forward(i,j) = (f_plus - f0) / h;
        
        % Central
        x_minus = x0_norm; 
        x_minus(i) = max(0, x0_norm(i) - h);
        f_minus = optim(x_minus, lb, ub, data, mda_options);
        all_dfdx_central(i,j) = (f_plus - f_minus) / (2*h);
    end
    
    % Plotting (using the new matrix indexing)
    subplot(2, 3, i);
    loglog(h_range, abs(all_dfdx_forward(i,:)), 'r--', 'LineWidth', 1.1); hold on;
    loglog(h_range, abs(all_dfdx_central(i,:)), 'b-', 'LineWidth', 1.5);
    grid on; title(['Sensitivity: ', var_names{i}]);
    xlabel('Step Size h'); ylabel(['|df/d', var_names{i}, '|']);
    hold on;

    xline(1e-3, 'k:', 'MDA Tol');
    xline(5e-3, 'g--', 'Chosen h');
end


