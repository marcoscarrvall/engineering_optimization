% Point 1 analysis
x_p1 = [0.02462, 0.004035];
h = 1e-8;

% Evaluate constraints to confirm activity
gx_p1 = springcon1(x_p1);
% Expected result: g(1), g(2), g(3) approx 0

% Gradients at Point 1
fx_p1 = springobj1(x_p1);
grad_f_p1 = [(springobj1([x_p1(1)+h, x_p1(2)]) - fx_p1)/h; ...
             (springobj1([x_p1(1), x_p1(2)+h]) - fx_p1)/h];

gx1_h = springcon1([x_p1(1)+h, x_p1(2)]);
gx2_h = springcon1([x_p1(1), x_p1(2)+h]);
grad_g_p1 = [(gx1_h - gx_p1)/h; (gx2_h - gx_p1)/h];

% Setup System for active constraints g1, g2, g3
% grad_g1*mu1 + grad_g2*mu2 + grad_g3*mu3 = -grad_f
A_p1 = [grad_g_p1(:,1), grad_g_p1(:,2), grad_g_p1(:,3)];
b_p1 = -grad_f_p1;

% Solve using pseudo-inverse due to dependency between g2 and g3
mu_p1 = pinv(A_p1) * b_p1;

% 4. Display Results
fprintf('Results for Intersection g1, g2, g3:\n');
fprintf("g1(x_p1) = %f\n", gx_p1(1));
fprintf("g2(x_p1) = %f\n", gx_p1(2));
fprintf("g3(x_p1) = %f\n", gx_p1(3));
fprintf("g4(x_p1) = %f\n", gx_p1(4));
fprintf("g5(x_p1) = %f\n", gx_p1(5));
fprintf('mu1: %f\n', mu_p1(1));
fprintf('mu2: %f\n', mu_p1(2));
fprintf('mu3: %f\n', mu_p1(3));
fprintf('Stationarity residual: %e\n', norm(A_p1 * mu_p1 + grad_f_p1));