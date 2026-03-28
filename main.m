x0 = [235, 5 ]; % [V, BPR]

coefficient = 0.1;

lb = [x0(1) * (1 - coefficient), x0(2) * (1 - coefficient)];
ub = [x0(1) * (1 + coefficient), x0(2) * (1 + coefficient)];

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...        
    'Display', 'iter-detailed', ... 
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6);


[x_opt, f_opt, exitflag] = fmincon(@(x) optim(x, TestAC_data), x0, [], [], [], [], lb, ub, @(x) constraints(x, TestAC_data), options);

fprintf('\n--- Optimization Results ---\n');
fprintf('Optimal Design Vector: [%s]\n', num2str(x_opt));
fprintf('Optimal Objective Value: %f\n', f_opt);

