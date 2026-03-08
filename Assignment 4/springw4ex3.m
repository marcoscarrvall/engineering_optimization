% Two variable valve spring problem - Exercise 4.3
% Solution of the unconstrained problem using the Matlab algorithm fminunc.

% 1. Initialization and Problem Visualization 
% (Re-using logic from springw4ex1.m to ensure Figure 1 is present)
clf, hold off, clear
format long

% Constant parameter values [cite: 18, 20, 21]
springparams1;
w = 1.0;
ktarget = 10000; 
frtarget = 300;

% Matrix of output values for combinations of design variables D and d [cite: 24]
D_range = [0.020:0.0005:0.040];
d_range = [0.002:0.00004:0.005];

for j = 1:length(d_range)
    for i = 1:length(D_range)
        % Analysis of valve spring for the contour plot
        [~,~,~,~,~,~,~,k,~,~,~,~,freq1] = ...
            springanalysis1(D_range(i), d_range(j), L0, L1, n, E, G, rho, Dv, h, p1, p2, nm, ncamfac, nne, matp, bldp);
        
        % Scaled objective function for visualization [cite: 14]
        fobj(j,i) = ((k-ktarget)/ktarget)^2 + w*((freq1-frtarget)/frtarget)^2; 
        stiffness(j,i) = k;
        freq(j,i) = freq1;
    end
end

% Create the Contour Plot [cite: 24, 25]
cc = [0.01 0.02 0.05];
contour(D_range, d_range, fobj, [cc 10*cc 100*cc 1000*cc 10000*cc 100000*cc 1000000*cc])
xlabel('Coil diameter D (m)'), ylabel('Wire diameter d (m)')
title('Figure 1: Valve Spring Optimization using fminunc (w = 1.0)')
hold on
grid on

% Plot the target lines (Intersection is the solution) [cite: 25]
contour(D_range, d_range, stiffness, [10000 10000], 'r', 'LineWidth', 1.5)
contour(D_range, d_range, freq, [300 300], 'g', 'LineWidth', 1.5)

% 2. Setup fminunc Optimization
% Initial design point [cite: 45]
xq = [0.022 0.0035];

% Configure options for medium-scale optimization and BFGS [cite: 63, 64]
% Note: 'quasi-newton' algorithm uses BFGS by default in fminunc.
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'HessUpdate', 'BFGS', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 20000, ... % Increased from 200
    'MaxIterations', 10000, ...          % Increased to allow more steps
    'TolFun', 1e-6); % Explicitly setting Tolfun as suggested [cite: 64]

% 3. Run the optimization
% Using an anonymous function to pass extra parameters to s_objw43 [cite: 60, 61]
[x_final, f_final, exitflag, output] = fminunc(...
    @(x) s_objw43(x, ktarget, frtarget, w), xq, options);

% 4. Plot the results [cite: 62, 63]
plot(xq(1), xq(2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Initial Point') % Initial point
plot(x_final(1), x_final(2), 'm*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Final Solution') % Final point
legend('Objective Contours', 'k = 10000 N/m', 'freq1 = 300 Hz', 'Initial Point', 'Final Solution')

% 5. Print Summary to Command Window
fprintf('\n--- Optimization Results ---\n');
fprintf('Initial Point: D = %.4f, d = %.4f\n', xq(1), xq(2));
fprintf('Final Point:   D = %.4f, d = %.4f\n', x_final(1), x_final(2));
fprintf('Final Objective Value: %e\n', f_final);
fprintf('Number of Iterations: %d\n', output.iterations);
fprintf('Number of Function Evaluations: %d\n', output.funcCount);