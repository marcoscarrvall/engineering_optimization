%% CLEARANCE.m
% Ground Clearance Constraint
%
% Computes the clearance between the bottom of the engine nacelle and the
% ground during static / take-off conditions, then checks it against the
% minimum allowable value.
%
% Geometry (all measured from ground level):
%
%   Ground ────────────────────────────────────── 0 m
%   Wing lower surface                           eng.wing_h_AGL
%   Fan centreline (= wing surface – pylon)      eng.wing_h_AGL – eng.pylon_h
%   Bottom of nacelle (= centreline – R_nac)     eng.wing_h_AGL – eng.pylon_h – R_nac
%   Clearance                                    above – 0
%
% Constraint is NORMALISED by clearance_min so the value is dimensionless
% and O(1) for use in gradient-based optimisers (e.g. fmincon):
%
%   c_cl = (clearance_min - clearance) / clearance_min  <=  0
%
% A value of +0.10 means clearance is 10% below the minimum required;
% a value of -0.20 means there is 20% more clearance than required.
%
% INPUTS:
%   state.D_nacelle   [m]   Nacelle outer diameter (from WATE)
%   eng.wing_h_AGL    [m]   Wing lower surface height above ground
%   eng.pylon_h       [m]   Pylon height (fan CL below wing lower surface)
%   con.clearance_min [m]   Minimum required clearance
%
% OUTPUTS:
%   c_cl      [-]   Normalised constraint  (<=0 satisfied, >0 violated)
%   violated  [bool]
%   margin    [-]   Normalised headroom  = -c_cl  (positive = safe)
%
% =========================================================================

function [c_cl] = clearance(state, eng, con)

    R_nacelle = state.D_nacelle / 2;

    % Vertical position of nacelle bottom above ground [m]
    clearance_val = eng.wing_h_AGL - eng.pylon_h - R_nacelle;

    % Store on state for post-processing / diagnostics
    state.clearance = clearance_val;   %#ok<NASGU>

    % Normalised inequality: g(x) = (clearance_min - clearance) / clearance_min <= 0
    c_cl     = (con.clearance_min - clearance_val) / con.clearance_min;
    violated = c_cl > 0;
    margin   = -c_cl;          % [-]  positive = safe headroom

    if violated
        fprintf('  STATUS        : *** VIOLATED (nacelle too large / mounted too low) ***\n');
    else
        fprintf('  STATUS        : Satisfied\n');
        fprintf('  Margin        : %+.6f [-]\n', margin);
    end

end