clear; clc; close all;

mda_options.tol = 1e-3;
mda_options.max_iter = 100;
mda_options.verbose = false;

data = A320data;

x = [220, 10.0, 1.55, 1.55, 22.0];

objective = optim(x, data, mda_options);

fprintf('Objective value: %.4f\n', objective);
fprintf("Range: %.2f km\n", -objective * data.ac.range / 1000);