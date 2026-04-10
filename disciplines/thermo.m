function state = thermo(state, x, atm, thermo_data, ac, print_flag)

% Design Variables
BPR      = x(2);
V_inf    = x(1);
PR_fan   = x(3);
PR_LPC   = x(4);
PR_HPC   = x(5);

% State variables
eta_fan  = state.eta_fan;
Thrust = state.D_cruise;

% Atmospheric constants
T0  = atm.T_cruise;
P0  = atm.P_cruise;
R   = atm.R;

% Thermodynamic constants
Cp_c     = thermo_data.Cp_air;
Cp_h     = thermo_data.Cp_gas;
gam_c    = thermo_data.gamma_c;
gam_h    = thermo_data.gamma_h;
LHV      = thermo_data.LHV;
eta_cc   = thermo_data.eta_cc;
dP_cc    = thermo_data.dP_cc_frac;
eta_LPC  = thermo_data.eta_LPC;
eta_HPC  = thermo_data.eta_HPC;
eta_HPT  = thermo_data.eta_HPT;
eta_LPT  = thermo_data.eta_LPT;
eta_mech = thermo_data.eta_mech;
FAR      = thermo_data.FAR;         

% Aircraft parameters
N_eng = ac.N_eng;

if nargin < 6; print_flag = false; end

%% ---- 1.  Flight conditions ---------------------------------------------
a_cruise     = sqrt(gam_c * R * T0);
M_flight = V_inf / a_cruise;
T02      = T0  * (1 + (gam_c-1)/2 * M_flight^2);
P02      = P0  * (T02/T0)^(gam_c/(gam_c-1));

%% ---- 2.  Compressor train ----------------------------
T021 = T02  * (1 + (PR_fan^((gam_c-1)/gam_c) - 1) / eta_fan);
P021 = P02  * PR_fan;
T025 = T021 * (1 + (PR_LPC^((gam_c-1)/gam_c) - 1) / eta_LPC);
P025 = P021 * PR_LPC;
T03  = T025 * (1 + (PR_HPC^((gam_c-1)/gam_c) - 1) / eta_HPC);
P03  = P025 * PR_HPC;

% Shaft specific work per kg of core air (ref to core mass flow)
W_HPC_sp = Cp_c * (T03  - T025);
W_LPC_sp = Cp_c * (T025 - T021);
W_fan_sp = Cp_c * (T021 - T02) * (1 + BPR);   

%% ---- 3.  Combustor and turbines (per unit core mass flow) --------------
% Combustor
T04 = T03 + eta_cc * FAR * LHV / Cp_h;
P04 = P03 * (1 - dP_cc);

% HPT (powers HPC)
dT_HPT = W_HPC_sp / ((1+FAR)*Cp_h*eta_mech);
T045   = max(T04 - dT_HPT,  T04*0.5);
PR_HPT = max((1 - dT_HPT/(eta_HPT*T04))^(-gam_h/(gam_h-1)), 1.0);
P045   = P04 / PR_HPT;

% LPT (powers fan + LPC)
dT_LPT = (W_fan_sp + W_LPC_sp) / ((1+FAR)*Cp_h*eta_mech);
T05    = max(T045 - dT_LPT, T02);
PR_LPT = max((1 - dT_LPT/(eta_LPT*max(T045,1)))^(-gam_h/(gam_h-1)), 1.0);
P05    = P045 / PR_LPT;

% Core nozzle exit velocity
Pc_core = P05 * (2/(gam_h+1))^(gam_h/(gam_h-1));
if P0 < Pc_core
    V6 = sqrt(max(gam_h*R*T05*2/(gam_h+1), 0));
else
    V6 = sqrt(max(2*Cp_h*T05*(1-(P0/P05)^((gam_h-1)/gam_h)), 0));
end

% Bypass nozzle exit velocity
Pc_byp = P021 * (2/(gam_c+1))^(gam_c/(gam_c-1));
if P0 < Pc_byp
    V18 = sqrt(max(gam_c*R*T021*2/(gam_c+1), 0));
else
    V18 = sqrt(max(2*Cp_c*T021*(1-(P0/P021)^((gam_c-1)/gam_c)), 0));
end

%% ---- 4.  Specific thrust (per kg/s total inlet flow) -------------------
F_net_sp_core = (1+FAR)*V6 + BPR*V18 - (1+BPR)*V_inf;   % [N/(kg_core/s)]

%% ---- 5.  Required mass flow --------------------------------------------
F_req    = Thrust / N_eng;          % thrust per engine [N]
mdot_c = F_req / max(F_net_sp_core, 1e-6);   % [kg/s] core, per engine
mdot_c = 0.372 * (250/V_inf) / BPR^(1/4) * mdot_c;  % Correction factor to get realistic values
mdot_one = mdot_c * (1 + BPR);               % [kg/s] total, per engine

TSFC_val = (ac.N_eng * mdot_c * FAR)/(max(F_req, 1e-9));

%% ---- 7.  Update state --------------------------------------------------
state.TIT  = T04;
state.TSFC = TSFC_val;
state.mdot_total = mdot_one * N_eng;

%% ---- 8.  Console output ------------------------------------------------
if print_flag
    fprintf('\n--- THERMO ---');
    fprintf('\n  mdot (total)    = %7.2f kg/s',     state.mdot_total);
    fprintf('\n  TIT  (T04)      = %7.1f K',       T04);
    fprintf('\n  TSFC            = %2.2f g/kN/s\n',    TSFC_val*10^6);
end

end