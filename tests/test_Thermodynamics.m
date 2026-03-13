% TEST_THERMODYNAMICS
%
% Regression test for Thermodynamics.m using parameters ported from the
% reference Python Params struct, extended with a realistic two-spool
% turbofan design point (BPR = 12.5).
%
% -------------------------------------------------------------------------
%  MAPPING FROM PYTHON PARAMS
%
%  Single efficiency → both spools:
%    compressor_eff = 0.87  →  eta_fan = eta_lpc = eta_hpc = 0.87
%    turbine_eff    = 0.90  →  eta_hpt = eta_lpt = 0.90
%
%  FAR back-calculated from delta_T = 900 K:
%    T3 is derived from the full fan → LPC → HPC compression chain.
%    T4 = T3 + 900 K is then enforced via the combustor energy balance:
%
%      FAR = (cp_gas*T4 - cp_air*T3) / (eta_b*LHV - cp_gas*T4)
%
%    so that the Thermodynamics solver reproduces exactly T4 = T3 + 900 K.
%
%  thrust_N is pre-computed so that the solver returns mdot_core = 16.3 kg/s
%  (equivalently mdot_total = 16.3 * (1 + BPR) = 220.05 kg/s).
%
% -------------------------------------------------------------------------
%  EXPECTED RESULTS
%    mdot_total  ≈ 220.0500 kg/s   → mdot_core = 16.3000 kg/s
%    TIT         ≈ 1645.51 K
%    TSFC        ≈ 18.29 g/kNs
% -------------------------------------------------------------------------

clear; clc;

%% --- data struct (ported from Python Params) ----------------------------

data.T_amb      = 216.5;          % K
data.p_amb      = 22632.0;        % Pa
data.R          = 287.0;          % J/kg/K
data.cp_air     = 1005.0;         % J/kg/K
data.kappa_air  = 1.4;
data.cp_gas     = 1150.0;         % J/kg/K
data.kappa_gas  = 1.33;
data.LHV        = 43e6;           % J/kg

% Mach 0.78 cruise → V0
data.V0 = 0.78 * sqrt(data.kappa_air * data.R * data.T_amb);  % 230.05 m/s

% Pressure ratios
data.inlet_pr   = 0.99;
data.comb_pr    = 0.96;

% Efficiencies — single value mapped to all spools
data.eta_fan    = 0.87;
data.eta_lpc    = 0.87;
data.eta_hpc    = 0.87;
data.eta_hpt    = 0.90;
data.eta_lpt    = 0.90;
data.eta_b      = 0.995;
data.eta_mech   = 0.995;
data.nozzle_eff = 0.98;

% Station 21 Mach number for A_core calculation
data.M21 = 0.5;

% --- FAR back-calculated to reproduce delta_T = 900 K -------------------
%
% Reproduce the compression chain to find T3, then solve the combustor
% energy balance for the FAR that gives T4 = T3 + 900 K.
%
%   cp_air*T3 + FAR*eta_b*LHV = (1 + FAR)*cp_gas*T4
%   => FAR = (cp_gas*T4 - cp_air*T3) / (eta_b*LHV - cp_gas*T4)

ka      = data.kappa_air;
Mach    = 0.78;
T0      = data.T_amb * (1 + (ka-1)/2 * Mach^2);
T2      = T0;
T21     = T2  * (1 + (1/data.eta_fan) * (1.35^((ka-1)/ka) - 1));  % fan_pr = 1.35
T25     = T21 * (1 + (1/data.eta_lpc) * (4.5^((ka-1)/ka)  - 1));  % lpc_pr = 4.5
T3      = T25 * (1 + (1/data.eta_hpc) * (5.5^((ka-1)/ka)  - 1));  % hpc_pr = 5.5
T4      = T3 + 900.0;

data.FAR = (data.cp_gas * T4 - data.cp_air * T3) / ...
           (data.eta_b  * data.LHV - data.cp_gas * T4);

fprintf('Compression chain:\n');
fprintf('  T2  = %8.4f K\n', T2);
fprintf('  T21 = %8.4f K\n', T21);
fprintf('  T25 = %8.4f K\n', T25);
fprintf('  T3  = %8.4f K\n', T3);
fprintf('  T4  = %8.4f K  (T3 + 900)\n', T4);
fprintf('  FAR = %.8f\n\n', data.FAR);

%% --- design vector ------------------------------------------------------

design_vec = [12.5, ...   % BPR
               1.35, ...  % fan_pr
               4.5,  ...  % lpc_pr
               5.5];      % hpc_pr

%% --- thrust target -------------------------------------------------------
%
% Pre-computed so that the solver returns mdot_core = 16.3 kg/s:
%   mdot_total = mdot_core * (1 + BPR) = 16.3 * 13.5 = 220.05 kg/s
%   thrust_N   = F_sp * mdot_total     = 113.2367 * 220.05 = 24917.73 N

thrust_N = 24917.7307;    % N

%% --- run solver ---------------------------------------------------------

[A_inlet, A_core, TIT, thrust, TSFC, mdot_total] = ...
    Thermodynamics(data, design_vec, thrust_N);

%% --- report results -----------------------------------------------------

BPR        = design_vec(1);
mdot_core  = mdot_total / (1 + BPR);

fprintf('--- Thermodynamics output ---\n');
fprintf('mdot_total  = %10.4f kg/s\n',   mdot_total);
fprintf('mdot_core   = %10.4f kg/s   (target: 16.3000)\n', mdot_core);
fprintf('TIT         = %10.4f K      (expected: ~1645.51)\n', TIT);
fprintf('thrust      = %10.4f N\n',   thrust);
fprintf('TSFC        = %10.4f g/kNs  (expected: ~18.29)\n', TSFC);
fprintf('A_inlet     = %10.4f m^2\n', A_inlet);
fprintf('A_core      = %10.4f m^2\n', A_core);

%% --- pass / fail --------------------------------------------------------

tol_mdot = 0.01;    % kg/s  on core flow
if abs(mdot_core - 16.3) <= tol_mdot
    fprintf('\nTEST PASSED: mdot_core = %.4f kg/s (within %.3f kg/s of 16.3)\n', ...
            mdot_core, tol_mdot);
else
    fprintf('\nTEST FAILED: mdot_core = %.4f kg/s (error = %.4f kg/s, tol = %.3f)\n', ...
            mdot_core, abs(mdot_core - 16.3), tol_mdot);
end