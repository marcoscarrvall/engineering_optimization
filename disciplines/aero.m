%% AERO.m
% Aerodynamics Module
%
% Computes the drag and lift changes that result from updated engine geometry
% and weight, then re-evaluates the aircraft lift-to-drag ratio.
%
% The drag model follows the reference formulation:
%   delta_CD0       — nacelle wetted-area increment relative to baseline
%   delta_CD_int    — interference drag scaling with (D_fan/D_ref)^2.5
%   delta_CD_induced— induced-drag change from the MTOW update
%
% INPUTS (via structs):
%   state.W_engine  [kg]   Current engine + pylon weight (one engine)
%   state.D_fan     [m]    Current fan tip diameter
%   state.D_nacelle [m]    Current nacelle outer diameter
%   state.L_engine  [m]    Current engine length
%   state.MTOW      [kg]   Current MTOW (will be updated here)
%   dv.V            [m/s]  Cruise true airspeed
%   ac              struct  Aircraft reference data
%   atm             struct  Atmospheric data
%   eng             struct  Engine reference data
%
% OUTPUTS (updated in state struct):
%   state.MTOW      [kg]   Updated MTOW
%   state.D_total   [N]    Total aircraft drag at cruise (both engines)
%   state.CL_CD     [-]    Updated cruise lift-to-drag ratio
%   state.CL        [-]    Current lift coefficient
% =========================================================================

function state = aero(state, ac, atm, dv, eng)

    V    = dv.V;
    rho  = atm.rho_cr;
    q    = 0.5 * rho * V^2;          % dynamic pressure [Pa]

    % =========================================================
    % 1.  UPDATED MTOW
    %     Engine weight change (×2 engines) propagates into OEW → MTOW
    % =========================================================
    delta_W_eng  = ac.N_engines * (state.W_engine - eng.W_engine_ref);
    MTOW_new     = ac.OEW + ac.payload + ac.fuel_mass + delta_W_eng;
    state.MTOW   = MTOW_new;

    % =========================================================
    % 2.  LIFT COEFFICIENT  (steady level flight: L = W)
    % =========================================================
    CL = (MTOW_new * atm.g) / (q * ac.S_ref);
    state.CL = CL;

    % =========================================================
    % 3.  PARASITIC DRAG INCREMENT  (nacelle wetted-area model)
    % =========================================================
    % Reference nacelle dimensions (CFM56-5B on A320)
    D_ref = eng.D_fan_ref;
    L_ref = eng.L_engine_ref;

    % Wetted-area ratio of current nacelle vs. reference
    S_wet_ratio   = (state.D_fan * state.L_engine) / (D_ref * L_ref);

    % Zero-lift drag increment from changed nacelle size
    delta_CD0     = 0.0020 * (S_wet_ratio - 1);

    % Interference drag — scales strongly with fan diameter growth
    delta_CD_int  = 0.0005 * (state.D_fan / D_ref)^2.5;

    delta_CD_parasitic = delta_CD0 + delta_CD_int;

    % =========================================================
    % 4.  INDUCED DRAG INCREMENT  (from MTOW change)
    % =========================================================
    % Change in CL relative to baseline
    CL_ref     = ac.CL_cr;
    delta_CL   = CL - CL_ref;

    % Total induced drag (using current CL, Oswald)
    CD_induced = CL^2 / (pi * ac.AR * ac.e);

    % =========================================================
    % 5.  TOTAL DRAG COEFFICIENT & DRAG FORCE
    % =========================================================
    CD_total   = ac.CD0_clean + delta_CD_parasitic + CD_induced;
    D_total    = CD_total * q * ac.S_ref;        % [N] both engines share this

    state.D_total = D_total;
    state.CL_CD   = CL / CD_total;

    fprintf('\n--- AERO ---\n');
    fprintf('  L/D             = %8.4f\n', state.CL_CD);
    fprintf('  Drag (total)    = %8.1f N\n', D_total);

end