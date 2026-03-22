%% run_MDA.m
% Sets up all input structs and calls the MDA converger.
% Edit the design variable dv.BPR (and optionally the sweep below) to
% explore different operating points.
% =========================================================================
clear; clc;

%% ── DESIGN VARIABLES ─────────────────────────────────────────────────────
% Only BPR is the active design variable in this MDA run.
% The remaining cycle DVs are held at their A320-class reference values.

dv.BPR    =  12.0;          % bypass ratio          [-]     (sweep this)
dv.PR_Fan =  1.60;         % fan pressure ratio    [-]
dv.PR_LPC =  1.80;         % LPC pressure ratio    [-]
dv.PR_HPC = 10.00;         % HPC pressure ratio    [-]
dv.V      = 230;           % cruise TAS            [m/s]   (≈ M0.78 @ FL350)

%% ── ATMOSPHERIC CONDITIONS  (ISA cruise, FL350) ─────────────────────────
atm.T_cr    = 218.8;       % static temperature    [K]
atm.P_cr    = 23842;       % static pressure       [Pa]
atm.rho_cr  = 0.3796;      % density               [kg/m³]
atm.gamma   = 1.40;        % ratio of specific heats (cold air)
atm.R       = 287.05;      % gas constant          [J/kg/K]
atm.g       = 9.80665;     % gravity               [m/s²]
atm.a_cr    = sqrt(atm.gamma * atm.R * atm.T_cr);   % speed of sound [m/s]

%% ── AIRCRAFT REFERENCE DATA ──────────────────────────────────────────────
ac.S_ref        = 122.4;   % wing reference area   [m²]       (A320-200)
ac.AR           =  9.39;   % aspect ratio          [-]
ac.e            =  0.82;   % Oswald efficiency     [-]
ac.CD0_clean    =  0.0230; % zero-lift drag coeff  [-]
ac.CL_cr        =  0.50;   % baseline cruise CL    [-]
ac.N_engines    =  2;      % number of engines     [-]
ac.OEW          = 42600;   % baseline OEW          [kg]   (without engines)
ac.payload      = 16600;   % design payload        [kg]
ac.fuel_mass    = 18800;   % block fuel            [kg]
ac.MTOW_ref     = 78000;   % reference MTOW        [kg]
ac.W_wing_ref   = 7700;    % reference wing weight [kg]

%% ── ENGINE REFERENCE DATA  (CFM56-5B) ───────────────────────────────────
eng.W_engine_ref  = 2800;  % reference engine+pylon weight [kg]  per engine
eng.D_fan_ref     = 1.735; % reference fan tip diameter    [m]
eng.L_engine_ref  = 3.50;  % reference engine length       [m]
eng.thrust_max    = 120e3; % sea-level static thrust (one) [N]
eng.rho_mat       = 4400;  % blade material density        [kg/m³]  (Ti alloy)

%% ── THERMODYNAMIC CONSTANTS ──────────────────────────────────────────────
thermo_data.Cp_air        = 1005;   % [J/kg/K]  cold stream
thermo_data.Cp_gas        = 1148;   % [J/kg/K]  hot stream
thermo_data.gamma_c       =  1.40;  % cold side
thermo_data.gamma_h       =  1.33;  % hot side
thermo_data.LHV           = 43.0e6; % lower heating value Jet-A  [J/kg]
thermo_data.eta_cc        =  0.995; % combustor efficiency
thermo_data.dP_cc_frac    =  0.04;  % combustor total-pressure loss fraction
thermo_data.eta_fan       =  0.88;
thermo_data.eta_LPC       =  0.87;
thermo_data.eta_HPC       =  0.86;
thermo_data.eta_HPT       =  0.89;
thermo_data.eta_LPT       =  0.91;
thermo_data.eta_mech      =  0.99;

%% ── WATE CALIBRATION (placeholder, kept for extensibility) ──────────────
wate_data = struct();

%% ── RUN SINGLE-POINT MDA ─────────────────────────────────────────────────
state = mda(dv, ac, atm, eng, thermo_data, wate_data, ...
            'conv_eps',   1e-4, ...
            'conv_niter', 40,   ...
            'display',    'iter+time');

%% ── PRINT CONVERGED SUMMARY ──────────────────────────────────────────────
fprintf('\n╔══════════════════════════════════════════╗\n');
fprintf('║       CONVERGED DESIGN POINT              ║\n');
fprintf('╠═══════════════════════════════════════════╣\n');
fprintf('║  BPR             = %5.2f                  ║\n', dv.BPR);
fprintf('║  OEW             = %7.0f kg               ║\n', state.mda.OEW);
fprintf('║  W_engine (one)  = %7.1f kg               ║\n', state.W_engine);
fprintf('║  mdot (both eng) = %7.2f kg/s             ║\n', state.mdot);
fprintf('║  L/D             = %7.4f                  ║\n', state.CL_CD);
fprintf('║  Drag            = %7.0f N                ║\n', state.D_total);
fprintf('╚═══════════════════════════════════════════╝\n\n');

%% ── OPTIONAL: BPR SWEEP ─────────────────────────────────────────────────
% Uncomment the block below to sweep BPR and plot key metrics.

% BPR_vec  = 4:0.5:12;
% results  = struct('BPR', num2cell(BPR_vec));
% 
% for i = 1:numel(BPR_vec)
%     dv_i      = dv;
%     dv_i.BPR  = BPR_vec(i);
%     s         = MDA(dv_i, ac, atm, eng, thermo_data, wate_data, ...
%                     'display', 'none');
%     results(i).TSFC      = s.TSFC;
%     results(i).W_engine  = s.W_engine;
%     results(i).L_D       = s.CL_CD;
%     results(i).OEW       = s.MDA.OEW;
%     results(i).D_fan     = s.D_fan;
% end
% 
% figure;
% tiledlayout(2,3);
% fields = {'TSFC','W_engine','L_D','OEW','D_fan'};
% ylabels = {'TSFC [kg/N/s]','W_{engine} [kg]','L/D [-]','OEW [kg]','D_{fan} [m]'};
% for k = 1:numel(fields)
%     nexttile;
%     plot(BPR_vec, [results.(fields{k})], 'o-', 'LineWidth', 1.5);
%     xlabel('BPR [-]'); ylabel(ylabels{k}); grid on;
%     title(ylabels{k});
% end
% sgtitle('BPR Sensitivity — MDA Converged Results');