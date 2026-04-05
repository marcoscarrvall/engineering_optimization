function state = thermo_simple(state, x, x_consts, atm, thermo_data, ac, print_flag)

BPR      = x(2);
V_inf    = x(1);

PR_fan   = x_consts.PR_fan;
PR_LPC   = x_consts.PR_LPC;
PR_HPC   = x_consts.PR_HPC;

Cp_c     = thermo_data.Cp_air;
Cp_h     = thermo_data.Cp_gas;
gam_c    = thermo_data.gamma_c;
gam_h    = thermo_data.gamma_h;
LHV      = thermo_data.LHV;
eta_cc   = thermo_data.eta_cc;
eta_fan  = thermo_data.eta_fan;
eta_LPC  = thermo_data.eta_LPC;
eta_HPC  = thermo_data.eta_HPC;
eta_HPT  = thermo_data.eta_HPT;
eta_LPT  = thermo_data.eta_LPT;
eta_mech = thermo_data.eta_mech;
FAR      = thermo_data.FAR;         

T0  = atm.T_cruise;
gam = atm.gamma;
R   = atm.R;


T = state.D_cruise;

if nargin < 7; print_flag = false; end

% Flow at BP nozzle (Assume non-choked)
PR_ram = 1 + (gam - 1)/2 * (V_inf / sqrt(gam*R*T0))^2;
T02 = T0 + V_inf^2 / (2*Cp_c);
tau_fan = 1 + 1/eta_fan * (PR_fan^((gam_c-1)/gam_c) - 1);
exp_nozzle = 1 - (1/ (PR_ram * PR_fan)^((gam_c-1)/gam_c));

V_bp = sqrt(2*Cp_c * T02 * tau_fan * exp_nozzle);

% Flow at core nozzle (Assume non-choked)

deltaT_fan = 1+ (PR_fan^((gam_c-1)/gam_c) - 1)/eta_fan;
deltaT_LPC = 1+ (PR_LPC^((gam_c-1)/gam_c) - 1)/eta_LPC;
deltaT_HPC = 1+ (PR_HPC^((gam_c-1)/gam_c) - 1)/eta_HPC;
deltaT_cc = eta_cc*FAR*LHV / Cp_h;

T04 = T02 * deltaT_fan * deltaT_LPC * deltaT_HPC + deltaT_cc;
W_HPT = Cp_c/ (Cp_h*eta_mech) * T02* deltaT_fan * deltaT_LPC * deltaT_HPC;
W_LPT = Cp_c / (Cp_h*eta_mech) * T02*deltaT_fan * (deltaT_LPC + (1+BPR));

core_energy_chain = T04 - W_HPT - W_LPT;

T021 = T02 * deltaT_fan;
T025 = T021 * deltaT_LPC;
T03  = T025 * deltaT_HPC;
T04 = T03 + (eta_cc * FAR * LHV) / Cp_h;
T045 = T04 - (Cp_c * (T03 - T025)) / (Cp_h * eta_mech);

work_fan = (1 + BPR) * (T021 - T02);
work_LPC = (T025 - T021);
T05 = T045 - (Cp_c * (work_fan + work_LPC)) / (Cp_h * eta_mech);

tau_g = gam_h / (gam_h - 1);
PR_HPT = (1 - (T04 - T045) / (T04 * eta_HPT)) ^ tau_g;
PR_LPT = (1 - (T045 - T05) / (T045 * eta_LPT)) ^ tau_g;

PR_core_exit = PR_ram * PR_fan * PR_LPC * PR_HPC * PR_HPT * PR_LPT;

V_core = sqrt(max(0, 2 * Cp_h * T05 * (1 - (1/PR_core_exit)^((gam_h-1)/gam_h))));


fprintf("T = %.2f N\n", T);
fprintf("V_inf = %.2f m/s\n", V_inf);
fprintf("V_bp = %.2f m/s\n", V_bp);
fprintf("V_core = %.2f m/s\n", V_core);
fprintf("BPR = %.2f\n", BPR);


mdot_c = T / (BPR * (V_bp - V_inf) + (V_core - V_inf));

mdot = mdot_c * (1 + BPR);

TSFC = FAR / ((V_core - V_inf) + BPR * (V_bp - V_inf));

%% ---- 7.  Update state --------------------------------------------------
state.TIT  = T04;
state.TSFC = TSFC;
state.mdot = mdot;

%% ---- 8.  Console output ------------------------------------------------
if print_flag
    fprintf('\n--- THERMO ---\n');
    fprintf('  TIT  (T04)      = %7.1f K\n',       T04);
    fprintf('  TSFC            = %.4e kg/N/s\n',    TSFC);
    fprintf('  mdot (total)    = %7.2f kg/s\n',     state.mdot);
    end

end