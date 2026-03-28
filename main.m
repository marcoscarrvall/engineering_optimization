x0 = [1.0, 1.0, 1.0, 1.0]; 

lb = [0.1, 0.1, 0.1, 0.1];
ub = [10.0, 10.0, 10.0, 10.0];

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...        
    'Display', 'iter-detailed', ... 
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6);


[x_opt, f_opt, exitflag] = fmincon(@objective_mdf, x0, [], [], [], [], lb, ub, @constraints, options);

fprintf('\n--- Optimization Results ---\n');
fprintf('Optimal Design Vector: [%s]\n', num2str(x_opt));
fprintf('Optimal Objective Value: %f\n', f_opt);