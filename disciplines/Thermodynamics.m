function [A_inlet, A_core, TIT, thrust, TSFC, mdot_total] = Thermodynamics(data, design_vec, thrust_N)
% THERMODYNAMICS  Two-spool two-stream turbofan cycle solver for MDA loops.
%
%   [A_inlet, A_core, TIT, thrust, TSFC, mdot_total] = Thermodynamics(data, design_vec, thrust_N)
%
%   Finds the total inlet mass flow that produces the required net thrust,
%   then computes the inlet capture area and core flow area from the
%   resulting thermodynamic state.
%
%   Architecture — two concentric spools:
%     LP spool : fan (full flow) + LPC (core only)  driven by LPT
%     HP spool : HPC (core only)                    driven by HPT
%
%   Combustor outlet temperature T4 is obtained from the fuel energy
%   balance using a specified fuel-to-air ratio (data.FAR):
%
%       T4 = [ mdot_core*cp_air*T3  +  mdot_fuel*eta_b*LHV ]
%            / [ mdot_gas * cp_gas ]
%
%   The cycle is linear in mass flow, so the solver runs at unit specific
%   mass flow (mdot_total = 1 kg/s) and scales to meet thrust_N.
%
% -------------------------------------------------------------------------
%   INPUTS
%     data        struct of constants.  Required fields:
%                   data.V0         freestream velocity           (m/s)
%                   data.T_amb      ambient static temperature    (K)
%                   data.p_amb      ambient static pressure       (Pa)
%                   data.R          specific gas constant, air    (J/kg/K)
%                   data.kappa_air  cp/cv ratio, air              (-)
%                   data.kappa_gas  cp/cv ratio, combustion gas   (-)
%                   data.cp_air     specific heat, air            (J/kg/K)
%                   data.cp_gas     specific heat, combustion gas (J/kg/K)
%                   data.inlet_pr   inlet total-pressure ratio    (-)
%                   data.eta_fan    fan isentropic efficiency     (-)
%                   data.eta_lpc    LPC isentropic efficiency     (-)
%                   data.eta_hpc    HPC isentropic efficiency     (-)
%                   data.eta_hpt    HPT isentropic efficiency     (-)
%                   data.eta_lpt    LPT isentropic efficiency     (-)
%                   data.eta_b      combustor efficiency          (-)
%                   data.eta_mech   mechanical (bearing) eff.     (-)
%                   data.comb_pr    combustor total-pressure ratio(-)
%                   data.nozzle_eff nozzle isentropic efficiency  (-)
%                   data.LHV        fuel lower heating value      (J/kg)
%                   data.FAR        fuel-to-air ratio             (-)
%
%     design_vec  [4x1] or [1x4]:
%                   design_vec(1) = bypass ratio              (BPR,    -)
%                   design_vec(2) = fan pressure ratio        (fan_pr, -)
%                   design_vec(3) = LPC pressure ratio        (lpc_pr, -)
%                   design_vec(4) = HPC pressure ratio        (hpc_pr, -)
%
%     thrust_N    required net thrust                          (N)
%
%   OUTPUTS
%     A_inlet     inlet capture area  = mdot_total / (rho0_static * V0) (m^2)
%                   rho0_static is the ambient STATIC density; this is the
%                   standard streamtube capture area used for nacelle sizing.
%     A_core      core flow area at fan exit / LPC face (station 21)    (m^2)
%                   A_core = mdot_core / (rho21 * V21)
%                   where rho21 and V21 are evaluated at the STATIC
%                   conditions assuming Mach = 0.5 at station 21.
%                   NOTE: station 21 Mach is assumed = data.M21 if provided,
%                   otherwise defaults to 0.5 (typical LPC face value).
%     TIT         turbine inlet temperature T4                (K)
%     thrust      net thrust  (equals thrust_N by construction)         (N)
%     TSFC        thrust-specific fuel consumption            (g/kNs)
%     mdot_total  total inlet mass flow that meets thrust_N   (kg/s)
%
% -------------------------------------------------------------------------
%
%   Station numbering (ARP755):
%     0    ambient
%     2    inlet exit                     (total)
%     21   fan exit / bypass split point  (total)
%     25   LPC exit                       (total)
%     3    HPC exit                       (total)
%     4    combustor exit / HPT inlet     (TIT = T4)
%     45   HPT exit / LPT inlet           (total)
%     5    LPT exit                       (total)
%     8    core nozzle throat/exit
%     18   bypass nozzle throat/exit
% -------------------------------------------------------------------------

    BPR    = design_vec(1);
    fan_pr = design_vec(2);
    lpc_pr = design_vec(3);
    hpc_pr = design_vec(4);
    V0     = data.V0;
    ka     = data.kappa_air;

    % --- ambient total conditions ----------------------------------------
    Mach = V0 / sqrt(ka * data.R * data.T_amb);
    T0   = data.T_amb * (1 + (ka-1)/2 * Mach^2);
    p0   = data.p_amb * (1 + (ka-1)/2 * Mach^2) ^ (ka/(ka-1));

    % --- inlet (station 2) -----------------------------------------------
    T2 = T0;
    p2 = data.inlet_pr * p0;

    % --- fan (station 2 -> 21) — full mass flow --------------------------
    T21 = T2 * (1 + (1/data.eta_fan) * (fan_pr^((ka-1)/ka) - 1));
    p21 = fan_pr * p2;

    % --- specific mass flows (normalised to mdot_total = 1 kg/s) ---------
    mdot_core_sp   = 1 / (1 + BPR);    % kg/s core   per kg/s total
    mdot_bypass_sp = BPR / (1 + BPR);  % kg/s bypass per kg/s total

    % --- bypass nozzle: fan exit (21) -> ambient (18) -------------------
    F_bypass_sp = nozzle_thrust(p21, T21, data.p_amb, ka, ...
                                data.nozzle_eff, data.cp_air, ...
                                mdot_bypass_sp, V0, data.R);

    % --- core stream -----------------------------------------------------
    [F_core_sp, TIT, fuel_flow_sp] = core_cycle( ...
        lpc_pr, hpc_pr, T2, T21, p21, mdot_core_sp, V0, data);

    % --- specific net thrust (N per kg/s of total inlet flow) ------------
    F_sp = F_core_sp + F_bypass_sp;

    if F_sp <= 0
        error('Thermodynamics:noThrust', ...
              'Cycle produces non-positive specific thrust (%.3f N*s/kg). ' ...
              'Check design_vec and data.', F_sp);
    end

    % --- scale to required thrust ----------------------------------------
    mdot_total = thrust_N / F_sp;
    mdot_core  = mdot_core_sp  * mdot_total;
    thrust     = thrust_N;                    % equals thrust_N by construction

    % --- TSFC ------------------------------------------------------------
    fuel_flow_total = fuel_flow_sp * mdot_total;        % kg/s actual
    TSFC            = fuel_flow_total / thrust_N * 1e6; % g/kNs

    % =====================================================================
    %   AREA OUTPUTS
    % =====================================================================

    % --- A_inlet : freestream capture area (station 0, static) ----------
    %   The inlet capture streamtube is defined by the freestream:
    %     A_inlet = mdot_total / (rho0_static * V0)
    %   rho0_static = p_amb / (R * T_amb)   [ambient STATIC density]

    rho0_static = data.p_amb / (data.R * data.T_amb);
    A_inlet     = mdot_total / (rho0_static * V0);

    % --- A_core : core annulus area at fan exit / LPC face (station 21) --
    %   Flow at station 21 is at total conditions (T21, p21).
    %   A Mach number is needed to convert to static conditions and velocity.
    %   Use data.M21 if provided; otherwise default to 0.5 (typical design).

    if isfield(data, 'M21')
        M21 = data.M21;
    else
        M21 = 0.5;
    end

    %   Static temperature and pressure at station 21:
    T21_static = T21 / (1 + (ka-1)/2 * M21^2);
    p21_static = p21 / (1 + (ka-1)/2 * M21^2) ^ (ka/(ka-1));

    %   Flow velocity at station 21:
    V21 = M21 * sqrt(ka * data.R * T21_static);

    %   Static density at station 21:
    rho21 = p21_static / (data.R * T21_static);

    %   Core annulus area (only the core fraction of the total fan face):
    A_core = mdot_core / (rho21 * V21);

end

function [F, T4, fuel_flow] = core_cycle( ...
        lpc_pr, hpc_pr, T2, T21, p21, mdot_core, V0, data)
%CORE_CYCLE  Two-spool core thermodynamic cycle.
%
%   Spool power balances (quantities at specific/unit total mass flow):
%     HP spool:  W_hpc              = W_hpt / eta_mech
%     LP spool:  W_fan_full + W_lpc = W_lpt / eta_mech
%
%   W_fan_full is the fan work on the TOTAL flow (core + bypass) because
%   the fan is mounted on the LP shaft.  Since mdot_total = 1 kg/s by
%   normalisation, W_fan_full = 1 * cp_air * (T21 - T2).

    ka = data.kappa_air;
    kg = data.kappa_gas;

    % -- LPC (station 21 -> 25) ------------------------------------------
    T25   = T21 * (1 + (1/data.eta_lpc) * (lpc_pr^((ka-1)/ka) - 1));
    p25   = lpc_pr * p21;
    W_lpc = mdot_core * data.cp_air * (T25 - T21);

    % -- HPC (station 25 -> 3) -------------------------------------------
    T3    = T25 * (1 + (1/data.eta_hpc) * (hpc_pr^((ka-1)/ka) - 1));
    p3    = hpc_pr * p25;
    W_hpc = mdot_core * data.cp_air * (T3 - T25);

    % -- combustor (station 3 -> 4): energy balance for T4 ---------------
    %
    %   FAR = mdot_fuel / mdot_core  (given via data.FAR)
    %
    %   Energy balance:
    %     mdot_fuel * eta_b * LHV = mdot_gas*cp_gas*T4 - mdot_core*cp_air*T3
    %
    %   Rearranged:
    %     T4 = (mdot_core*cp_air*T3 + mdot_fuel*eta_b*LHV) / (mdot_gas*cp_gas)

    FAR       = data.FAR;
    fuel_flow = mdot_core * FAR;          % specific fuel mass flow
    mdot_gas  = mdot_core + fuel_flow;    % = mdot_core * (1 + FAR)

    T4 = (mdot_core * data.cp_air * T3 + fuel_flow * data.eta_b * data.LHV) ...
         / (mdot_gas * data.cp_gas);
    p4 = data.comb_pr * p3;

    % -- HP turbine (station 4 -> 45): drives HPC on HP spool ------------
    W_hpt = W_hpc / data.eta_mech;           % work removed from gas
    T45   = T4 - W_hpt / (mdot_gas * data.cp_gas);

    exp_hpt = 1 - (1/data.eta_hpt) * (1 - T45/T4);
    if exp_hpt <= 0
        F = -1e9;  fuel_flow = 0;  return
    end
    p45 = exp_hpt ^ (kg/(kg-1)) * p4;

    % -- LP turbine (station 45 -> 5): drives fan (full flow) + LPC ------
    %
    %   mdot_total_sp = 1 kg/s by normalisation convention.
    %   Fan work on full specific flow:
    %     W_fan_full = 1 [kg/s] * cp_air * (T21 - T2)
    %   LP spool balance:
    %     W_lpt = (W_fan_full + W_lpc) / eta_mech

    W_fan_full = 1.0 * data.cp_air * (T21 - T2);   % mdot_total_sp = 1 kg/s
    W_lpt      = (W_fan_full + W_lpc) / data.eta_mech;

    T5 = T45 - W_lpt / (mdot_gas * data.cp_gas);

    exp_lpt = 1 - (1/data.eta_lpt) * (1 - T5/T45);
    if exp_lpt <= 0
        F = -1e9;  fuel_flow = 0;  return
    end
    p5 = exp_lpt ^ (kg/(kg-1)) * p45;

    % -- core nozzle (station 5 -> 8) ------------------------------------
    F = nozzle_thrust(p5, T5, data.p_amb, kg, ...
                      data.nozzle_eff, data.cp_gas, mdot_gas, V0, data.R);
end

function F = nozzle_thrust(pa, Ta, p_amb, kappa, eff, cp, mdot, V0, R)
%NOZZLE_THRUST  Net thrust from a convergent nozzle (choked or unchoked).
%
%   pa, Ta   total pressure & temperature at nozzle entry
%   p_amb    ambient static pressure
%   kappa    ratio of specific heats
%   eff      nozzle isentropic efficiency
%   cp       specific heat of working gas
%   mdot     mass flow through nozzle
%   V0       aircraft velocity (ram-drag reference)
%   R        specific gas constant

    if pa <= p_amb
        F = mdot * (0 - V0);
        return
    end

    % critical pressure ratio (nozzle chokes when pa/p_amb >= cpr)
    cpr = (1 - (1/eff) * (kappa-1)/(kappa+1)) ^ (-kappa/(kappa-1));

    if (pa / p_amb) < cpr
        % --- unchoked: fully expanded to p_amb ---------------------------
        Tb = Ta * (1 + eff * ((p_amb/pa)^((kappa-1)/kappa) - 1));
        dT = Ta - Tb;
        if dT <= 0
            F = mdot * (0 - V0);
            return
        end
        Vb = sqrt(2 * cp * dT);
        F  = mdot * (Vb - V0);
    else
        % --- choked: exit plane at critical conditions -------------------
        pb  = pa / cpr;
        Tb  = Ta * 2 / (kappa + 1);
        Vb  = sqrt(kappa * R * Tb);
        rho = pb / (R * Tb);
        Ab  = mdot / (rho * Vb);
        F   = mdot * (Vb - V0) + (pb - p_amb) * Ab;
    end
end