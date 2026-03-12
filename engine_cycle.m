function [comp_pr, TIT, TSFC] = engine_cycle(data, design_vec, thrust_N)
% ENGINE_CYCLE  Two-stream turbofan cycle solver for MDA loops.
%
%   [comp_pr, TIT, TSFC] = engine_cycle(data, design_vec, thrust_N)
%
%   Finds the compressor pressure ratio that produces the required thrust,
%   given a bypass ratio and freestream velocity as design variables.
%   The combustor temperature rise (delta_T) and fan pressure ratio
%   (fan_pr) are held constant in data.
%
% -------------------------------------------------------------------------
%   INPUTS
%     data        struct of constants — load with data.m
%                 Must include data.fan_pr (fan pressure ratio)
%     design_vec  [2x1] or [1x2] vector:
%                   design_vec(1) = bypass ratio          (BPR, -)
%                   design_vec(2) = freestream velocity   (V0,  m/s)
%     thrust_N    required net thrust (N)
%
%   OUTPUTS
%     comp_pr     compressor pressure ratio that meets thrust_N  (-)
%     TIT         turbine inlet temperature T4                   (K)
%     TSFC        thrust-specific fuel consumption               (g/kNs)
% -------------------------------------------------------------------------
%
%   Station numbering:
%     0   ambient
%     2   inlet exit          (total)
%     21  fan exit            (total) — bypass splits here
%     3   compressor exit     (total)
%     4   combustor exit / turbine inlet  (TIT = T4)
%     5   turbine exit        (total)
%     8   core nozzle exit
%     18  bypass nozzle exit
% -------------------------------------------------------------------------

    BPR = design_vec(1);
    V0  = design_vec(2);

    % --- ambient total conditions ----------------------------------------
    Mach = V0 / sqrt(data.kappa_air * data.R * data.T_amb);
    T0   = data.T_amb * (1 + (data.kappa_air - 1)/2 * Mach^2);
    p0   = data.p_amb * (1 + (data.kappa_air - 1)/2 * Mach^2) ...
               ^ (data.kappa_air / (data.kappa_air - 1));

    % --- inlet (station 2) -----------------------------------------------
    T2 = T0;
    p2 = data.inlet_pr * p0;

    % --- fan (station 21) — applied to full mass flow --------------------
    T21 = T2 * (1 + (1/data.eta_fan) * ...
                (data.fan_pr^((data.kappa_air-1)/data.kappa_air) - 1));
    p21 = data.fan_pr * p2;

    % --- mass flow split after fan ---------------------------------------
    mdot_core   = data.mdot_total / (1 + BPR);
    mdot_bypass = data.mdot_total - mdot_core;

    % --- bypass stream: expands from fan exit (station 21 → 18) ---------
    F_bypass = nozzle_thrust(p21, T21, data.p_amb, data.kappa_air, ...
                             data.nozzle_eff, data.cp_air, ...
                             mdot_bypass, V0, data.R);

    % --- residual for fzero ----------------------------------------------
    function res = residual(comp_pr)
        F_core = core_thrust(comp_pr, T21, p21, mdot_core, V0, data);
        res    = (F_core + F_bypass) - thrust_N;
    end

    % --- bracket ---------------------------------------------------------
    pr_lo = 1.5;
    pr_hi = 40.0;

    F_lo = core_thrust(pr_lo, T21, p21, mdot_core, V0, data) + F_bypass;
    F_hi = core_thrust(pr_hi, T21, p21, mdot_core, V0, data) + F_bypass;

    if (F_lo - thrust_N) * (F_hi - thrust_N) > 0
        error('engine_cycle: target thrust %.1f N is outside achievable range [%.1f, %.1f] N.', ...
              thrust_N, min(F_lo, F_hi), max(F_lo, F_hi));
    end

    % --- solve -----------------------------------------------------------
    opts    = optimset('TolX', 1e-8, 'Display', 'off');
    comp_pr = fzero(@residual, [pr_lo, pr_hi], opts);

    % --- final state -----------------------------------------------------
    [~, TIT, fuel_flow, F_core] = core_thrust(comp_pr, T21, p21, mdot_core, V0, data);
    F_total = F_core + F_bypass;
    TSFC    = fuel_flow / F_total * 1e6;    % g/kNs

end


% =========================================================================
%   LOCAL HELPER FUNCTIONS
% =========================================================================

function [F, T4, fuel_flow, F_out] = core_thrust(comp_pr, T21, p21, mdot_core, V0, data)
% Core stream: compressor → combustor → turbine → nozzle.
% Inlet is fan exit (T21, p21). Compressor PR is on top of fan PR.

    % -- compressor (station 21 → 3) --------------------------------------
    T3     = T21 * (1 + (1/data.eta_c) * ...
                    (comp_pr^((data.kappa_air-1)/data.kappa_air) - 1));
    p3     = comp_pr * p21;
    W_comp = mdot_core * data.cp_air * (T3 - T21);

    % -- combustor (fixed delta_T) ----------------------------------------
    T4        = T3 + data.delta_T;
    p4        = data.comb_pr * p3;
    fuel_flow = mdot_core * data.cp_gas * data.delta_T / (data.eta_b * data.LHV);
    mdot_gas  = mdot_core + fuel_flow;

    % -- turbine: extracts W_comp / eta_mech ------------------------------
    W_turb  = W_comp / data.eta_mech;
    T5      = T4 - W_turb / (mdot_gas * data.cp_gas);
    exp_arg = 1 - (1/data.eta_t) * (1 - T5/T4);
    if exp_arg <= 0
        F     = -1e9;
        F_out = F;
        return
    end
    p5 = exp_arg ^ (data.kappa_gas / (data.kappa_gas - 1)) * p4;

    % -- core nozzle (station 5 → 8) --------------------------------------
    F     = nozzle_thrust(p5, T5, data.p_amb, data.kappa_gas, ...
                          data.nozzle_eff, data.cp_gas, mdot_gas, V0, data.R);
    F_out = F;
end


function F = nozzle_thrust(pa, Ta, p_amb, kappa, eff, cp, mdot, V0, R)
% Net thrust from a convergent nozzle — choked or unchoked.

    if pa <= p_amb
        F = mdot * (0 - V0);
        return
    end

    cpr = (1 - (1/eff) * (kappa-1)/(kappa+1)) ^ (-kappa/(kappa-1));

    if (pa / p_amb) < cpr
        % unchoked
        Tb = Ta * (1 + eff * ((p_amb/pa)^((kappa-1)/kappa) - 1));
        dT = Ta - Tb;
        if dT <= 0
            F = mdot * (0 - V0);
            return
        end
        Vb = sqrt(2 * cp * dT);
        F  = mdot * (Vb - V0);
    else
        % choked
        pb  = pa / cpr;
        Tb  = Ta * (2 / (kappa + 1));
        Vb  = sqrt(kappa * R * Tb);
        rho = pb / (R * Tb);
        Ab  = mdot / (rho * Vb);
        F   = mdot * (Vb - V0) + (pb - p_amb) * Ab;
    end
end