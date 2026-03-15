%% TIT.m
% Turbine Inlet Temperature Constraint
%
% Ensures TIT does not exceed the material design limit.
%
% Constraint is NORMALISED by TIT_max so the value is dimensionless
% and O(1) for use in gradient-based optimisers (e.g. fmincon):
%
%   c_TIT = (TIT - TIT_max) / TIT_max  <=  0
%
% A value of +0.01 means TIT is 1% above the limit;
% a value of -0.05 means TIT has 5% headroom.
%
% INPUTS:
%   state.TIT     [K]  Turbine inlet temperature (from THERMO)
%   con.TIT_max   [K]  Maximum allowable TIT (from Data.m)
%
% OUTPUTS:
%   c_TIT     [-]   Normalised constraint  (<=0 satisfied, >0 violated)
%   violated  [bool]
%   margin    [-]   Normalised headroom  = -c_TIT  (positive = safe)
%
% =========================================================================

function [c_TIT, violated, margin] = tit(state, con)

    TIT     = state.TIT;
    TIT_max = con.TIT_max;

    % Normalised inequality: g(x) = (TIT - TIT_max) / TIT_max <= 0
    c_TIT   = (TIT - TIT_max) / TIT_max;

    violated = c_TIT > 0;
    margin   = -c_TIT;          % [-]  positive = safe headroom

    fprintf('\n--- TIT Constraint ---\n');
    fprintf('  TIT         = %.1f K\n',  TIT);
    fprintf('  TIT_max     = %.1f K\n',  TIT_max);
    fprintf('  c_TIT       = %+.6f [-]  (normalised)\n', c_TIT);
    fprintf('  Margin      = %+.6f [-]  (= %.1f K)\n', margin, margin * TIT_max);
    if violated
        fprintf('  STATUS      : *** VIOLATED ***\n');
    else
        fprintf('  STATUS      : Satisfied\n');
    end

end