clear; clc; close all;

lb = [200, 5.0];
ub = [240, 14.0];

start_points = [   
    240, 14.0;  
    240, 5.1; 
    200, 5.0;   
    200, 14.0;
];

var_names        = {'V (m/s)', 'BPR'};
constraint_names = {'Clearance', 'TIT Limit', 'Tip Mach'};

mda_options.tol      = 1e-3;
mda_options.max_iter = 100;
mda_options.verbose  = false;
data = A320data;

max_slp_iter = 60;
M_penalty    = 1e6;
h_step       = 5e-3;
conv_tol     = 1e-4;

figure('Color', 'w', 'Name', 'Multi-Start Optimization Paths', 'Units', 'normalized');
hold on; grid on;

n_grid = 25; 
n_grid = 25; 
v_grid_axis = linspace(lb(1), ub(1), n_grid);
bpr_grid_axis = linspace(lb(2), ub(2), n_grid);
[V_mesh, BPR_mesh] = meshgrid(v_grid_axis, bpr_grid_axis);

% Pre-allocate matrices for speed
Z_obj = zeros(n_grid, n_grid);
G1 = zeros(n_grid, n_grid);
G2 = zeros(n_grid, n_grid);
G3 = zeros(n_grid, n_grid);

fprintf('Generating background topology...\n');
for r = 1:n_grid        % Loop through rows
    for c = 1:n_grid    % Loop through columns
        % Get the physical values from the mesh at this specific grid point
        V_val = V_mesh(r, c);
        BPR_val = BPR_mesh(r, c);
        
        % Normalize the current point for the functions
        x_norm = normalize_vars([V_val, BPR_val], lb, ub);
        
        % Calculate Objective
        Z_obj(r, c) = simplified_optim(x_norm, lb, ub, data, mda_options);
        
        % Calculate Constraints
        g_tmp = simplified_constraints(x_norm, lb, ub, data, mda_options);
        G1(r, c) = g_tmp(1); 
        G2(r, c) = g_tmp(2); 
        G3(r, c) = g_tmp(3);
    end
end

% Plot Background
[C, h_cont] = contour(V_mesh, BPR_mesh, Z_obj, 20, 'k', 'LineWidth', 0.5);
clabel(C, h_cont, 'FontSize', 8);
[~, h1] = contour(V_mesh, BPR_mesh, G1, [0 0], 'g', 'LineWidth', 1.5);
[~, h2] = contour(V_mesh, BPR_mesh, G2, [0 0], 'b', 'LineWidth', 1.5);
[~, h3] = contour(V_mesh, BPR_mesh, G3, [0 0], 'c', 'LineWidth', 1.5);


colors = lines(size(start_points, 1)); 

final_objectives = zeros(size(start_points, 1), 1);

for s = 1:size(start_points, 1)
    x_curr_phys = start_points(s, :);
    x_curr = normalize_vars(x_curr_phys, lb, ub);
    move_limit = 0.1;
    path_history = x_curr_phys; 
        
    for k = 1:max_slp_iter
        f_curr = simplified_optim(x_curr, lb, ub, data, mda_options);
        [g_curr, ~] = simplified_constraints(x_curr, lb, ub, data, mda_options);
        
        [grad_f, grad_g] = get_gradients(x_curr, lb, ub, data, mda_options, h_step);
        
        % Solve Linear Program
        f_lp = [grad_f(:); M_penalty];
        A_lp = [grad_g, -ones(size(grad_g,1), 1)];
        b_lp = -g_curr;
        lb_dx = max(-move_limit, 0 - x_curr);
        ub_dx = min(move_limit, 1 - x_curr);
        
        opt_lp = optimoptions('linprog', 'Display', 'none');
        [res, ~, exitflag] = linprog(f_lp, A_lp, b_lp, [], [], [lb_dx(:); 0], [ub_dx(:); inf], opt_lp);
        
        if exitflag ~= 1, break; end
        
        dx = res(1:2)';
        x_next = x_curr + dx;
        
        % Acceptance Logic
        f_next = simplified_optim(x_next, lb, ub, data, mda_options);
        [g_next, ~] = simplified_constraints(x_next, lb, ub, data, mda_options);
        
        if max(g_next) > 1e-6 && max(g_next) > max(g_curr) || (max(g_curr) <= 1e-6 && f_next > f_curr)
            move_limit = move_limit * 0.5;
        else
            x_curr = x_next;
            path_history = [path_history; denormalize_vars(x_curr, lb, ub)];
            if norm(dx) < conv_tol, break; end
            if norm(dx) >= 0.99 * move_limit, move_limit = min(0.2, move_limit * 1.2); end
        end
    end
    final_objectives(s) = simplified_optim(x_curr, lb, ub, data, mda_options);
    
    plot(path_history(:,1), path_history(:,2), '-o', 'Color', colors(s,:), ...
         'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', colors(s,:));
    plot(path_history(1,1), path_history(1,2), 'o', 'MarkerEdgeColor', colors(s,:), 'MarkerSize', 4, 'LineWidth', 1);
end

% Final Plot Formatting
xlabel('Velocity (V)'); ylabel('Bypass Ratio (BPR)');
title('SLP Multi-Start Convergence Analysis');
legend([h1, h2, h3], {'Clearance Limit', 'TIT Limit', 'Mach Limit'}, 'Location', 'northeastoutside');
axis([lb(1) ub(1) lb(2) ub(2)]);

% Give final objective values for each starting point
fprintf('Final objective values for each starting point:\n');
for s = 1:size(start_points, 1)
    fprintf('Starting Point %d: %.4f\n', s, final_objectives(s));
end

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
    n = length(x);
    [g0, ~] = simplified_constraints(x, lb, ub, data, opts);
    grad_g = zeros(length(g0), n);
    grad_f = zeros(n, 1);
    for i = 1:n
        xp = x; xp(i) = min(1, x(i)+h);
        xm = x; xm(i) = max(0, x(i)-h);
        fp = simplified_optim(xp, lb, ub, data, opts);
        fm = simplified_optim(xm, lb, ub, data, opts);
        grad_f(i) = (fp - fm) / (xp(i) - xm(i));
        gp = simplified_constraints(xp, lb, ub, data, opts);
        gm = simplified_constraints(xm, lb, ub, data, opts);
        grad_g(:,i) = (gp - gm) / (xp(i) - xm(i));
    end
end