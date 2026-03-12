function data = data()
% DATA  Central data file for the MDA loop.
%
%   data = data()
%
%   All constants and fixed parameters for the analysis live here.
%   Group them by discipline. The engine_cycle function reads from
%   the fields defined in the [ENGINE] section — do not rename those.

% =========================================================================
%   MISSION / FLIGHT CONDITIONS
% =========================================================================
data.altitude_m     = 11000;        % cruise altitude                   (m)
data.T_amb          = 216.5;        % ambient static temperature        (K)
data.p_amb          = 22632.0;      % ambient static pressure           (Pa)
data.rho_amb        = data.p_amb / (287.0 * data.T_amb);  % air density (kg/m^3)

% =========================================================================
%   AIRCRAFT / WEIGHTS
% =========================================================================
data.g              = 9.80665;      % gravitational acceleration        (m/s^2)
data.W_start        = 600000;       % aircraft weight at start of cruise (N)  (~71 t)
% data.MTOW           = ;           % max take-off weight               (kg)
% data.OEW            = ;           % operating empty weight            (kg)
% data.payload        = ;           % design payload                    (kg)
% data.W_fuel         = ;           % fuel weight                       (kg)

% =========================================================================
%   AERODYNAMICS
% =========================================================================
data.CL_CD          = 18.0;         % cruise lift-to-drag ratio         (-)
% data.S_ref          = ;           % reference wing area               (m^2)
% data.AR             = ;           % aspect ratio                      (-)
% data.e              = ;           % Oswald efficiency factor          (-)
% data.CD0            = ;           % zero-lift drag coefficient        (-)
% data.CL_cruise      = ;           % cruise lift coefficient           (-)

% =========================================================================
%   ENGINE — required by engine_cycle.m (do not rename these fields)
% =========================================================================
data.mdot_total     = 103.0;        % total inlet mass flow             (kg/s)
data.inlet_pr       = 0.99;         % inlet pressure recovery           (-)
data.fan_pr         = 1.6;          % fan pressure ratio                (-)
data.eta_fan        = 0.92;         % fan isentropic efficiency         (-)
data.comb_pr        = 0.96;         % combustor pressure ratio          (-)
data.delta_T        = 900.0;        % combustor temperature rise T4-T3  (K)

data.eta_c          = 0.87;         % compressor isentropic efficiency  (-)
data.eta_t          = 0.90;         % turbine isentropic efficiency     (-)
data.eta_mech       = 0.995;        % mechanical shaft efficiency       (-)
data.eta_b          = 0.995;        % combustor efficiency              (-)
data.nozzle_eff     = 0.98;         % nozzle efficiency                 (-)

data.R              = 287.0;        % specific gas constant — air       (J/kg/K)
data.cp_air         = 1005.0;       % specific heat — cold section      (J/kg/K)
data.kappa_air      = 1.40;         % ratio of specific heats — air     (-)
data.cp_gas         = 1150.0;       % specific heat — hot section       (J/kg/K)
data.kappa_gas      = 1.33;         % ratio of specific heats — gas     (-)
data.LHV            = 43.0e6;       % fuel lower heating value          (J/kg)

% =========================================================================
%   STRUCTURES
% =========================================================================
% data.n_ult          = ;           % ultimate load factor              (-)
% data.t_over_c       = ;           % wing thickness-to-chord ratio     (-)
% data.material_rho   = ;           % structural material density       (kg/m^3)
% data.sigma_allow    = ;           % allowable stress                  (Pa)

% =========================================================================
%   PERFORMANCE
% =========================================================================
data.W_fuel         = 120000;       % fixed fuel weight                 (N)  (~12.2 t)
% data.range          = ;           % design range (use breguet_fuel.m) (m)
data.V_cruise       = 230;          % cruise true airspeed              (m/s) (~Mach 0.78)
% data.Mach_cruise    = ;           % cruise Mach number                (-)
% data.SFC_target     = ;           % target SFC                        (g/kNs)

end