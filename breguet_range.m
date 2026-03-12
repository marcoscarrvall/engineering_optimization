function range = breguet_range(data, design_vec, W_start, thrust_N)
% BREGUET_RANGE  Range achievable for a fixed fuel weight.
%
%   range = breguet_range(data, design_vec, W_start, thrust_N)
%
%   Breguet range equation solved for R:
%
%       R = (V0 / (g * TSFC_SI)) * (L/D) * ln(W_start / W_end)
%
%   where  W_end = W_start - W_fuel  and TSFC_SI [kg/(N*s)] = TSFC [g/kNs] * 1e-6
%
% -------------------------------------------------------------------------
%   INPUTS
%     data        struct from data.m — must include:
%                   data.W_fuel    fixed fuel weight           (N)
%                   data.CL_CD     lift-to-drag ratio          (-)
%                   data.g         gravitational acceleration  (m/s^2)
%     design_vec  [BPR, V0]  — passed directly to engine_cycle
%     W_start     aircraft weight at start of cruise          (N)
%     thrust_N    cruise thrust per engine                    (N)
%
%   OUTPUT
%     range       achievable range                            (m)
% -------------------------------------------------------------------------
 
    [~, ~, TSFC_g_kNs] = engine_cycle(data, design_vec, thrust_N);
 
    TSFC_SI = TSFC_g_kNs * 1e-6;
 
    V0    = design_vec(2);
    W_end = W_start - data.W_fuel;
 
    range = (V0 / (data.g * TSFC_SI)) * data.CL_CD * log(W_start / W_end);
 
end