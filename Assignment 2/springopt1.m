clear; clc; close all;

springparams1;

freq1lb = ncamfac * (nm / 2) / 60;

% Design variable range
D_range = 0.02:0.0005:0.04; % Mean winding diameter (m)
d_range = 0.002:0.0001:0.005; % Wire diameter (m)

% For plotting later on
OBJ = zeros(length(d_range), length(D_range));
G1 = OBJ; G2 = OBJ; G3 = OBJ; G4 = OBJ; G5 = OBJ; G6 = OBJ;

for j = 1:length(d_range)
    for i = 1:length(D_range)
        [svol,smass,bvol,matc,manc,Lmin,L2,k,F1,F2,Tau1,Tau2,freq1] = ...
            springanalysis1(D_range(i), d_range(j), L0, L1, n, E, G, rho, Dv, h, p1, p2, nm, ncamfac, nne, matp, bldp);
        
        Av = (pi/4)*Dv^2;
        F1min = Av * p1;
        F2min = Av * p2;
        Tau12max = 600e6;

        OBJ(j,i) = smass;

        G1(j,i) = (Lmin - L2) / Lmin;            
        G2(j,i) = (F1min - F1) / F1min;         
        G3(j,i) = (F2min - F2) / F2min;          
        G4(j,i) = (Tau2 - Tau12max) / Tau12max;  
        G5(j,i) = (freq1lb - freq1) / freq1lb;   
    end
end

% Visualization of objective function
hold on
[C, h_obj] = contour(D_range, d_range, OBJ, 'k:');
clabel(C, h_obj);

% Visualization of constraints
contour(D_range, d_range, G1, [0.0 0.0], 'r', 'LineWidth', 2); % Lmin
contour(D_range, d_range, G2, [0.0 0.0], 'g', 'LineWidth', 2); % F1min
contour(D_range, d_range, G3, [0.0 0.0], 'b', 'LineWidth', 2); % F2min
contour(D_range, d_range, G4, [0.0 0.0], 'm', 'LineWidth', 2); % Tau2max
contour(D_range, d_range, G5, [0.0 0.0], 'c', 'LineWidth', 2); % freq1lb

xlabel('Mean winding diameter D [m]');
ylabel('Wire diameter d [m]');
title('Optimization Problem: Minimize smass');
legend('smass (Obj)', 'g1 (Lmin)', 'g2 (F1min)', 'g3 (F2min)', 'g4 (Tau2)', 'g5 (freq1)', ...
       'Location', 'northeastoutside');
grid on