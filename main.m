clear; clc; close all;

constraint_names = {'clearance', 'noise', 'TIT'};

% Initial guess and bounds
x0 = [235, 5]; 
x_consts.PR_fan = 1.7;
x_consts.PR_LPC = 2.6;
x_consts.PR_HPC = 6.1;


options_mda.tol = 1e-6;
options_mda.max_iter = 100;
options_mda.verbose = false;

coefficient = 4;
lb = x0 * (1 - coefficient);
ub = x0 * (1 + coefficient);

% Global storage
global optHistory current_c
optHistory.fval = []; optHistory.constr = []; optHistory.iter = [];

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'none', ... 
    'OutputFcn', @outfun, ... 
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6);

fprintf('\n%-10s %-12s %-12s %-12s\n', 'Iter', 'V (x1)', 'BPR (x2)', 'f(x)');
fprintf('-----------------------------------------------------------\n');

[x_opt, f_opt] = fmincon(@(x) optim(x, x_consts, TestAC_data, options_mda), x0, ...
                         [], [], [], [], lb, ub, ...
                         @(x) constraints_catcher(x, x_consts, TestAC_data, options_mda), options);

% --- PLOTTING SECTION ---
figure('Color', 'w', 'Name', 'Optimization Convergence');

% Top Plot: Objective
subplot(2,1,1);
plot(optHistory.iter, optHistory.fval, '-bo', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on; ylabel('Objective f(x)'); title('Minimize Objective');

% Bottom Plot: Constraints
subplot(2,1,2);
if ~isempty(optHistory.constr)
    h = plot(optHistory.iter, optHistory.constr, '-s', 'LineWidth', 1.2);
    hold on;
    yline(0, 'r--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
    grid on; ylabel('Value (c <= 0)'); xlabel('Iteration');
    title('Constraint Satisfaction');
    
    % --- Apply the Labels ---
    % If you have more constraints than labels, it defaults to 'C1', 'C2', etc.
    num_found = size(optHistory.constr, 2);
    if length(constraint_names) >= num_found
        legend(h, constraint_names(1:num_found), 'Location', 'bestoutside');
    else
        legend(h, 'Location', 'bestoutside'); 
    end
end

% --- INTERNAL FUNCTIONS ---

function [c, ceq] = constraints_catcher(x, x_consts, data, opts)
    global current_c
    [c, ceq] = constraints(x, x_consts, data, opts);
    current_c = c; % Capture current values for the logger
end

function stop = outfun(x, optimValues, state)
    stop = false;
    global optHistory current_c
    if strcmp(state, 'iter')
        fprintf('%-10d %-12.4f %-12.4f %-12.6f\n', ...
            optimValues.iteration, x(1), x(2), optimValues.fval);
        
        optHistory.iter = [optHistory.iter; optimValues.iteration];
        optHistory.fval = [optHistory.fval; optimValues.fval];
        optHistory.constr = [optHistory.constr; current_c(:)'];
    end
end