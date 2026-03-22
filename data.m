%% Data.m
% Stores all constants and initial computed parameters for the
% aero-engine MDO optimizer.
% Reference aircraft  :  Airbus A320  /  CFM56-5B turbofan
% =========================================================================

%% =========================================================
%  1. ATMOSPHERE  (ISA at cruise  —  35 000 ft / 10 668 m)
% =========================================================
atm.T_sl   = 288.15;          % [K]     Sea-level temperature
atm.P_sl   = 101325;          % [Pa]    Sea-level pressure
atm.rho_sl = 1.225;           % [kg/m3] Sea-level density
atm.gamma  = 1.4;             % [-]     Ratio of specific heats (air)
atm.R      = 287.05;          % [J/kg/K]Gas constant (air)
atm.g      = 9.80665;         % [m/s2]  Gravitational acceleration

atm.h_cr   = 10668;           % [m]     35 000 ft cruise altitude
atm.T_cr   = 218.81;          % [K]     ISA cruise static temperature
atm.P_cr   = 23842;           % [Pa]    ISA cruise static pressure
atm.rho_cr = 0.3796;          % [kg/m3] ISA cruise density
atm.a_cr   = sqrt(atm.gamma * atm.R * atm.T_cr);  % [m/s] speed of sound

%% =========================================================
%  2. DESIGN VECTOR  —  initial / reference values
% =========================================================
dv.V      = 230;              % [m/s]   Cruise TAS  (~M 0.78 at 35 kft)
dv.BPR    = 6.0;              % [-]     Bypass ratio  (CFM56-5B ref)
dv.PR_HPC = 9.5;              % [-]     HPC pressure ratio
dv.PR_LPC = 2.6;              % [-]     LPC (booster) pressure ratio
dv.PR_Fan = 1.7;              % [-]     Fan pressure ratio
dv.OPR    = dv.PR_Fan * dv.PR_LPC * dv.PR_HPC;  % [-] Overall PR

%% =========================================================
%  3. AIRCRAFT GEOMETRY & WEIGHTS  (A320 reference)
% =========================================================
ac.MTOW       = 77000;        % [kg]   Max Take-Off Weight
ac.OEW        = 42600;        % [kg]   Operating Empty Weight
ac.fuel_mass  = 20100;        % [kg]   Max usable fuel
ac.payload    = 14300;        % [kg]   Design payload (150 pax × ~95 kg)
ac.W_wing_ref = 7700;         % [kg]   Reference wing weight (for scaling)

ac.S_ref      = 122.4;        % [m2]   Wing reference area
ac.b          = 34.1;         % [m]    Wing span
ac.AR         = ac.b^2 / ac.S_ref;   % [-]  Aspect ratio ~ 9.5
ac.e          = 0.82;         % [-]    Oswald efficiency factor
ac.sweep_LE   = 27.1;         % [deg]  Wing LE sweep

% Aerodynamic reference (cruise baseline)
ac.CL_cr      = 0.56;         % [-]    Baseline cruise CL
ac.CD0_clean  = 0.0220;       % [-]    Zero-lift parasitic drag (clean)
ac.CL_CD_ref  = 17.5;         % [-]    Reference cruise L/D
ac.N_engines  = 2;            % [-]    Number of engines

%% =========================================================
%  4. ENGINE REFERENCE DATA  (CFM56-5B, per engine)
% =========================================================
eng.N_engines    = 2;
eng.mdot_ref     = 350;       % [kg/s]  Total mass flow reference (both engines)
eng.T_sl_ref     = 120000;    % [N]     Sea-level static thrust per engine
eng.thrust_max   = 120000;    % [N]     Used in pylon weight calculation
eng.TSFC_ref     = 1.65e-5;   % [kg/N/s]Cruise TSFC reference
eng.W_engine_ref = 2350;      % [kg]    Reference engine + pylon weight (one engine)
eng.D_fan_ref    = 1.735;     % [m]     Reference fan diameter (CFM56-5B)
eng.L_engine_ref = 3.50;      % [m]     Reference engine length
eng.rho_mat      = 4430;      % [kg/m3] Titanium alloy (fan blade material)

% Geometric installation parameters
eng.pylon_h      = 0.40;      % [m]     Pylon height  (fan CL to wing lower surface)
eng.wing_h_AGL   = 1.80;      % [m]     Wing lower surface height above ground (static)

%% =========================================================
%  5. THERMODYNAMIC CONSTANTS
% =========================================================
therm.Cp_air     = 1005;     % [J/kg/K] Cold-stream specific heat
therm.Cp_gas     = 1150;     % [J/kg/K] Hot-gas specific heat
therm.gamma_c    = 1.40;     % [-]      Cold stream gamma
therm.gamma_h    = 1.33;     % [-]      Hot stream gamma
therm.LHV        = 43.2e6;   % [J/kg]   Lower heating value Jet-A
therm.eta_cc     = 0.995;    % [-]      Combustion efficiency
therm.dP_cc_frac = 0.05;     % [-]      Fractional total pressure loss in combustor (5 %)
therm.eta_fan    = 0.92;     % [-]      Fan isentropic efficiency
therm.eta_LPC    = 0.90;     % [-]      LPC isentropic efficiency
therm.eta_HPC    = 0.87;     % [-]      HPC isentropic efficiency
therm.eta_HPT    = 0.91;     % [-]      HPT isentropic efficiency
therm.eta_LPT    = 0.93;     % [-]      LPT isentropic efficiency
therm.eta_mech   = 0.99;     % [-]      Mechanical transmission efficiency

% Fixed fuel-to-air ratio (per kg core air)
% CFM56-5B cruise FAR ~ 0.0270  (gives TIT ~ 1480 K at cruise conditions)
therm.FAR        = 0.0270;   % [-]      Fuel-to-air ratio (constant, design point)

%% =========================================================
%  6. CONSTRAINT THRESHOLDS
% =========================================================
con.TIT_max       = 1700;     % [K]   Max TIT — single-crystal Ni limit
con.clearance_min = 0.45;     % [m]   Min ground clearance (fan bottom to ground)
con.M_tip_max     = 1.50;     % [-]   Max fan-tip relative Mach (noise l

%% =========================================================
%  8. MISSION PROFILE
% =========================================================
mis.R_design  = 3300e3;       % [m]    Design range (3 300 km, typical A320)
mis.h_cr      = atm.h_cr;     % [m]    Cruise altitude
mis.reserve_f = 0.05;         % [-]    Reserve fuel fraction

%% =========================================================
%  9. INITIAL MDA STATE  (updated iteratively)
% =========================================================
state.mdot      = eng.mdot_ref;     % [kg/s]  Initial total mass flow (both engines)
state.W_engine  = eng.W_engine_ref; % [kg]    Initial engine + pylon weight (one engine)
state.MTOW      = ac.MTOW;          % [kg]    Initial MTOW
state.D_total   = 60000 * 2;        % [N]     Initial drag estimate (60 kN/engine × 2)
state.CL_CD     = ac.CL_CD_ref;     % [-]     Initial L/D
state.TSFC      = eng.TSFC_ref;     % [kg/N/s]Initial TSFC
state.D_fan     = eng.D_fan_ref;    % [m]     Initial fan diameter
state.D_nacelle = eng.D_fan_ref * 1.10;  % [m] Initial nacelle diameter
state.L_engine  = eng.L_engine_ref; % [m]     Initial engine length
state.A_fan     = 0;                % [m2]    (filled by WATE)
state.N_stages  = 0;                % [-]     (filled by WATE)
state.TIT       = 0;                % [K]     (filled by THERMO)
state.TSFC      = 0;                % [kg/N/s](filled by THERMO)
state.range     = 0;                % [m]     (filled by BREGUET)
state.M_tip     = 0;                % [-]     (filled by NOISE)
state.CL        = 0;                % [-]     (filled by AERO)
state.clearance = 0;                % [m]     (filled by CLEARANCE)

fprintf('Data.m loaded — A320 / CFM56-5B reference.\n');
fprintf('  Design range : %.0f km\n', mis.R_design/1e3);
fprintf('  MTOW         : %.0f kg\n', ac.MTOW);
fprintf('  BPR / OPR    : %.1f / %.1f\n', dv.BPR, dv.OPR);
fprintf('  Fixed FAR    : %.4f\n', therm.FAR);