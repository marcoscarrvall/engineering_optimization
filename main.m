clear; clc; close all;

% Initial Vector [V, BPR, PR_fan, PR_LPC, PR_HPC]
x0 = [230, 8.0, 1.55, 1.55, 22.0];
lb = [200, 5.0, 1.1, 1.1, 10.0];
ub = [235, 14.0, 1.6, 1.8, 25.0];

% No constraints for now
A = []; b = []; Aeq = []; beq = [];

% Optimization options
optimizer_options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'none');


% MDA options
mda_options.tol = 1e-3;
mda_options.max_iter = 100;
mda_options.verbose = false;

% Load data
data = A320data;

% Run optimization
[x, fval] = fmincon(@(x) optim(x, data, mda_options), x0, A, b, Aeq, beq, lb, ub, [], optimizer_options);