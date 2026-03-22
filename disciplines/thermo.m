function state = thermo(state, dv, atm, thermo_data, verbose)
%% THERMO.m  —  Thermodynamic cycle analysis (design-point turbofan)
%
% DROP-IN REPLACEMENT for the original fixed-FAR thermo.m.
% Same 4-argument signature (verbose defaults to true).
% Call as thermo(state, dv, atm, thermo_data, false) to suppress output.
%
% Key change: FAR is iterated via bisection so that the net thrust per
% engine equals state.D_total / 2.  TIT is therefore a genuine function
% of BPR and operating point.  state.mdot is also updated so WATE can
% resize the engine on the next MDA iteration.
%
% Coupling equations solved here:
%   F_net(FAR, mdot, BPR) = state.D_total / 2          [thrust match]
%   state.mdot = (state.D_total / Fsp) * 1             [mdot from Fsp]
%   (mdot is total for both engines; Fsp = F_net/mdot_one_engine)
%
% Station numbering (ARP755 / SAE):
%   2→021 fan   025 LPC exit   3 HPC exit
%   4 combustor (TIT)   45 HPT exit   5 LPT exit
%   6 core nozzle exit   18 bypass nozzle exit
% =========================================================================

%% ---- 0.  Unpack --------------------------------------------------------
BPR      = dv.BPR;
PR_fan   = dv.PR_Fan;
PR_LPC   = dv.PR_LPC;
PR_HPC   = dv.PR_HPC;
V_inf    = dv.V;

Cp_c     = thermo_data.Cp_air;
Cp_h     = thermo_data.Cp_gas;
gam_c    = thermo_data.gamma_c;
gam_h    = thermo_data.gamma_h;
LHV      = thermo_data.LHV;
eta_cc   = thermo_data.eta_cc;
dP_cc    = thermo_data.dP_cc_frac;
eta_fan  = thermo_data.eta_fan;
eta_LPC  = thermo_data.eta_LPC;
eta_HPC  = thermo_data.eta_HPC;
eta_HPT  = thermo_data.eta_HPT;
eta_LPT  = thermo_data.eta_LPT;
eta_mech = thermo_data.eta_mech;

T0  = atm.T_cr;
P0  = atm.P_cr;
gam = atm.gamma;
R   = atm.R;

N_eng = 2;
if nargin < 5; verbose = true; end

%% ---- 1.  Flight conditions ---------------------------------------------
a_cr     = sqrt(gam * R * T0);
M_flight = V_inf / a_cr;
T02      = T0  * (1 + (gam_c-1)/2 * M_flight^2);
P02      = P0  * (T02/T0)^(gam_c/(gam_c-1));

%% ---- 2.  Compressor train (FAR-independent) ----------------------------
T021 = T02  * (1 + (PR_fan^((gam_c-1)/gam_c) - 1) / eta_fan);
P021 = P02  * PR_fan;
T025 = T021 * (1 + (PR_LPC^((gam_c-1)/gam_c) - 1) / eta_LPC);
P025 = P021 * PR_LPC;
T03  = T025 * (1 + (PR_HPC^((gam_c-1)/gam_c) - 1) / eta_HPC);
P03  = P025 * PR_HPC;

% Shaft specific work per kg of CORE air
W_HPC_sp = Cp_c * (T03  - T025);               % high spool compressor
W_LPC_sp = Cp_c * (T025 - T021);               % low  spool compressor
W_fan_sp = Cp_c * (T021 - T02) * (1 + BPR);    % fan (all flow, ref to core mass)

%% ---- 3.  Required thrust per engine ------------------------------------
F_req = state.D_total / N_eng;   % [N]

%% ---- 4.  Cycle function (nested, closes over compressor state) ---------
    function [F_net, V6, V18, T04, T045, T05, Fsp, TSFC_v, mdot_one_v] = ...
            cycle(FAR, mdot_one)
        % mdot_one : total mass flow through one engine [kg/s]
        mdot_c = mdot_one / (1 + BPR);   % core stream
        mdot_b = mdot_one * BPR/(1+BPR); % bypass stream

        % Combustor
        T04 = (Cp_c*T03 + FAR*LHV*eta_cc) / (Cp_h*(1+FAR));
        P04 = P03 * (1 - dP_cc);

        % HPT (high spool: powers HPC)
        dT_HPT = W_HPC_sp / ((1+FAR)*Cp_h*eta_mech);
        T045   = max(T04 - dT_HPT,  T04*0.5);
        PR_HPT = max((1 - dT_HPT/(eta_HPT*T04))^(-gam_h/(gam_h-1)), 1.0);
        P045   = P04 / PR_HPT;

        % LPT (low spool: powers fan + LPC)
        dT_LPT = (W_fan_sp + W_LPC_sp) / ((1+FAR)*Cp_h*eta_mech);
        T05    = max(T045 - dT_LPT, T02);
        PR_LPT = max((1 - dT_LPT/(eta_LPT*max(T045,1)))^(-gam_h/(gam_h-1)), 1.0);
        P05    = P045 / PR_LPT;

        % Core nozzle
        Pc_core = P05 * (2/(gam_h+1))^(gam_h/(gam_h-1));
        if P0 < Pc_core
            V6 = sqrt(max(gam_h*R*T05*2/(gam_h+1), 0));
        else
            V6 = sqrt(max(2*Cp_h*T05*(1-(P0/P05)^((gam_h-1)/gam_h)), 0));
        end

        % Bypass nozzle
        Pc_byp = P021 * (2/(gam_c+1))^(gam_c/(gam_c-1));
        if P0 < Pc_byp
            V18 = sqrt(max(gam_c*R*T021*2/(gam_c+1), 0));
        else
            V18 = sqrt(max(2*Cp_c*T021*(1-(P0/P021)^((gam_c-1)/gam_c)), 0));
        end

        % Net thrust and metrics
        F_net = mdot_c*(1+FAR)*V6 - mdot_c*V_inf + mdot_b*V18 - mdot_b*V_inf;
        Fsp   = F_net / max(mdot_one, 1e-9);
        TSFC_v      = (mdot_c*FAR) / max(F_net, 1e-9);
        mdot_one_v  = mdot_one;
    end

%% ---- 5.  Solve for mdot and FAR simultaneously -------------------------
%
% The two unknowns are coupled:
%   (a) FAR: for a given mdot, bisect to find FAR such that F_net = F_req
%   (b) mdot: after finding FAR, update mdot = F_req / Fsp
%
% This is a 2×2 fixed-point iteration (outer: mdot, inner: FAR bisection).
% It converges quickly because F_net is nearly linear in both variables.

mdot_one = state.mdot / N_eng;   % start from current state value

for outer = 1:20   % outer loop: converge mdot
    mdot_one_old = mdot_one;

    % ---- Inner bisection on FAR for current mdot ------------------------
    FAR_lo = 0.0;
    FAR_hi = 0.10;
    [F_lo,~,~,~,~,~,~,~,~] = cycle(FAR_lo, mdot_one);
    [F_hi,~,~,~,~,~,~,~,~] = cycle(FAR_hi, mdot_one);

    if F_lo >= F_req
        % Fan alone overshoots — mdot too large, will be corrected below
        FAR_sol = FAR_lo;
    elseif F_hi < F_req
        % Can't reach thrust even at max FAR — use ceiling
        FAR_sol = FAR_hi;
    else
        for k = 1:60
            FAR_mid = 0.5*(FAR_lo + FAR_hi);
            [F_mid,~,~,~,~,~,~,~,~] = cycle(FAR_mid, mdot_one);
            if abs(F_mid - F_req) < 0.1; break; end
            if F_mid < F_req; FAR_lo = FAR_mid; else; FAR_hi = FAR_mid; end
        end
        FAR_sol = FAR_mid;
    end

    % ---- Update mdot from Fsp -------------------------------------------
    [~,~,~,~,~,~, Fsp_cur,~,~] = cycle(FAR_sol, mdot_one);
    mdot_one_new = F_req / max(Fsp_cur, 1e-3);

    % Relaxed update to aid convergence
    mdot_one = 0.5*mdot_one_new + 0.5*mdot_one;

    if abs(mdot_one - mdot_one_old)/max(mdot_one_old,1) < 1e-4
        break;
    end
end

%% ---- 6.  Final evaluation at converged (FAR, mdot) --------------------
[~, V6, V18, T04, T045, T05, Fsp, TSFC_val, ~] = cycle(FAR_sol, mdot_one);

state.TIT  = T04;
state.TSFC = TSFC_val;
state.mdot = mdot_one * N_eng;   % update coupling variable for WATE

%% ---- 7.  Console output ------------------------------------------------
if verbose
    fprintf('\n--- THERMO ---\n');
    fprintf('  TIT  (T04)      = %7.1f K\n',    T04);
    fprintf('  TSFC            = %.4e kg/N/s\n',    TSFC_val);
    fprintf('  mdot (total)    = %7.2f kg/s\n',     state.mdot);
end

end