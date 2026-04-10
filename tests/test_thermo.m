clear; clc;

%% 1. Design variables
V_inf    = 240;     % Velocity in m/s (~Mach 0.78 at 35k ft)
BPR      = 15.0;     % Bypass Ratio
PR_fan   = 1.6;     % Fan Pressure Ratio
PR_LPC   = 1.8;     % Low Pressure Compressor PR
PR_HPC   = 25.0;    % High Pressure Compressor PR

x = [V_inf, BPR, PR_fan, PR_LPC, PR_HPC];

data = A320data; 

options.tol = 1e-3;
options.max_iter = 100;
options.verbose = true;


V_range = linspace(200, 235, 20);
BPR_range = linspace(5, 14, 20);
for i = 1:length(V_range)
    for j = 1:length(BPR_range)
        x(1) = V_range(i);
        x(2) = BPR_range(j);
        state = mda(x, data, options);
        objective(i,j) = state.objective;
        constraint1(i,j) = state.clearance_constraint;
        constraint2(i,j) = state.TIT_constraint;
        constraint3(i,j) = state.M_tip_constraint;
    end
end

figure;

[C, h] = contour(V_range, BPR_range, real(objective'), 20, 'k');
clabel(C, h); 
grid on;

xlabel('Velocity V (m/s)');
ylabel('Bypass Ratio (BPR)');
title('Contour Plot of Objective Function');

[min_val, idx] = min(real(objective(:)));
[v_idx, bpr_idx] = ind2sub(size(objective), idx);
hold on;

plot(V_range(v_idx), BPR_range(bpr_idx), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
legend('Objective Contours', 'Optimum Design');

hold on;
[C1, h1] = contour(V_range, BPR_range, real(constraint1'), [0 0], 'r', 'LineWidth', 0.5);
clabel(C1, h1, 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');

hold on;
[C2, h2] = contour(V_range, BPR_range, real(constraint2'), [0 0], 'b', 'LineWidth', 0.5);
clabel(C2, h2, 'FontSize', 10, 'Color', 'b', 'FontWeight', 'bold');

hold on; 
[C3, h3] = contour(V_range, BPR_range, real(constraint3'), [0 0], 'g', 'LineWidth', 0.5);
clabel(C3, h3, 'FontSize', 10, 'Color', 'g', 'FontWeight', 'bold');