classdef TestAC_data
% TESTAC_DATA  All fixed constants for the two-spool turbofan MDA test case.
%
%   Covers three groups:
%     1. Thermodynamic constants  (inputs to Thermodynamics.m)
%     2. Wing / aircraft constants (inputs to WATE.m)
%     3. Reference engine constants (baseline for WATE delta calculations)
%
%   Thermodynamics outputs (mdot_total, TIT, thrust, TSFC, fuel_flow,
%   A_inlet, A_core) are NOT stored here -- they are written into the
%   live data struct at runtime by the MDA loop before calling WATE.
%
%   Design variables (BPR, fan_pr, lpc_pr, hpc_pr) and thrust_N are also
%   NOT stored here -- they are supplied separately by the caller.
%
%   Usage:
%     data = TestAC_data.data_struct();

    properties (Constant)

        name = "testAC"

        % =================================================================
        %  1. THERMODYNAMIC CONSTANTS
        % =================================================================

        % --- Ambient / flight conditions ---------------------------------
        T_amb   = 216.5        % [K]    ISA tropopause static temperature
        p_amb   = 22632.0      % [Pa]   ISA tropopause static pressure
        Mach_cr = 0.78         % [-]    cruise Mach number
        V0      = 230.053327   % [m/s]  freestream velocity

        % --- Gas properties ----------------------------------------------
        R          = 287.0     % [J/kg/K]  specific gas constant, air
        cp_air     = 1005.0    % [J/kg/K]  specific heat, air
        kappa_air  = 1.4       % [-]       cp/cv, air
        cp_gas     = 1150.0    % [J/kg/K]  specific heat, combustion gas
        kappa_gas  = 1.33      % [-]       cp/cv, combustion gas
        LHV        = 43e6      % [J/kg]    fuel lower heating value

        % --- Pressure ratios ---------------------------------------------
        inlet_pr  = 0.99       % [-]  inlet total-pressure recovery
        comb_pr   = 0.96       % [-]  combustor total-pressure ratio

        % --- Component efficiencies --------------------------------------
        eta_fan    = 0.87      % [-]  fan isentropic efficiency
        eta_lpc    = 0.87      % [-]  LPC isentropic efficiency
        eta_hpc    = 0.87      % [-]  HPC isentropic efficiency
        eta_hpt    = 0.90      % [-]  HPT isentropic efficiency
        eta_lpt    = 0.90      % [-]  LPT isentropic efficiency
        eta_b      = 0.995     % [-]  combustor efficiency
        eta_mech   = 0.995     % [-]  mechanical efficiency
        nozzle_eff = 0.98      % [-]  nozzle isentropic efficiency

        % --- Combustor ---------------------------------------------------
        %   FAR pre-computed from delta_T = 900 K:
        %     T3 = 745.51 K  (fan_pr=1.35, lpc_pr=4.5, hpc_pr=5.5)
        %     T4 = 1645.51 K
        %     FAR = (cp_gas*T4 - cp_air*T3) / (eta_b*LHV - cp_gas*T4)
        FAR = 0.0279536484     % [-]  fuel-to-air ratio

        % --- Area calculation --------------------------------------------
        M21 = 0.5              % [-]  Mach number at LPC face (station 21)

        % =================================================================
        %  2. AIRCRAFT / WING CONSTANTS
        % =================================================================

        % --- Wing geometry -----------------------------------------------
        S_ref  = 122.6         % [m^2]  wing reference area
        AR     = 9.5           % [-]    wing aspect ratio
        e      = 0.85          % [-]    Oswald span efficiency
        sweep  = 25.0          % [deg]  wing quarter-chord sweep (LE sweep ~ 30 deg)
        t_c    = 0.12          % [-]    wing thickness-to-chord ratio (root mean)
        taper  = 0.25          % [-]    wing taper ratio lambda
        eta_eng = 0.35         % [-]    non-dim engine spanwise station (2*y/b)

        % --- Baseline aerodynamics (clean wing, no engine effect) --------
        CD0_ref  = 0.0200      % [-]  baseline zero-lift drag coefficient
        CL_ref   = 0.50        % [-]  reference cruise lift coefficient

        % --- Aircraft weights --------------------------------------------
        MTOW     = 77000       % [kg]  max take-off weight
        OEW_ref  = 42000       % [kg]  reference OEW (baseline engine)
        W_fuel   = 20000       % [kg]  fuel load at start of cruise

        % --- Mission parameters ------------------------------------------
        R_range      = 5500e3      % [m]   design range (5500 km)
        fuel_reserve = 0.05        % [-]   reserve fuel fraction (5%)

        % --- Engine count ------------------------------------------------
        n_eng    = 2           % [-]  number of engines

        % --- Material density for engine structure -----------------------
        rho_mat  = 4430        % [kg/m^3]  titanium alloy (fan / compressor discs)

        % =================================================================
        %  3. REFERENCE ENGINE  (CFM56-5B baseline for delta calculations)
        %     Used to compute delta_CD, delta_CL, delta_OEW relative to a
        %     known installation.  All values are per engine.
        % =================================================================

        % --- Reference geometry ------------------------------------------
        D_fan_ref  = 1.73      % [m]   fan tip diameter
        L_eng_ref  = 3.50      % [m]   overall engine length
        S_wet_ref  = 19.0      % [m^2] nacelle wetted area (= pi*D*L approx)

        % --- Reference weight (installed, per engine) --------------------
        W_eng_ref  = 2400      % [kg]  engine + pylon, per engine

        % --- Reference thrust (sea-level static, per engine) -------------
        thrust_max_ref = 120000  % [N]  max take-off thrust per engine

        % --- Reference OPR -----------------------------------------------
        OPR_ref = 27.0         % [-]  overall pressure ratio of reference engine

        % =================================================================
        %  4. CONSTRAINT THRESHOLDS
        % =================================================================

        % --- TIT (turbine inlet temperature) -----------------------------
        TIT_max    = 1800.0    % [K]   max allowable TIT (material limit)

        % --- Noise / fan tip Mach number ---------------------------------
        M_tip_max  = 1.05      % [-]   max fan tip relative Mach (noise/efficiency)

        % --- Ground clearance --------------------------------------------
        h_engine   = 2.10      % [m]   engine centreline height above ground
                               %       (pylon + wing dihedral, typical A320-class)
        clearance_min = 0.30   % [m]   minimum fan tip-to-ground clearance [m]
        tip_gap_frac  = 0.005  % [-]   min tip-to-nacelle gap as fraction of D_fan

        % --- Endurance ---------------------------------------------------
        t_max      = 14.0 * 3600   % [s]   max flight time (crew duty limit, 14 h)

    end

    methods (Static)

        function data = data_struct()
        % DATA_STRUCT  Return a plain struct for use by Thermodynamics.m and WATE.m.
        %
        %   Fields from Thermodynamics.m (mdot_total, TIT, thrust, TSFC,
        %   fuel_flow, A_inlet, A_core) are initialised to NaN here and
        %   must be filled by the MDA loop before calling WATE.

            % -- thermodynamic constants --
            data.T_amb      = TestAC_data.T_amb;
            data.p_amb      = TestAC_data.p_amb;
            data.Mach_cr    = TestAC_data.Mach_cr;
            data.V0         = TestAC_data.V0;
            data.R          = TestAC_data.R;
            data.cp_air     = TestAC_data.cp_air;
            data.kappa_air  = TestAC_data.kappa_air;
            data.cp_gas     = TestAC_data.cp_gas;
            data.kappa_gas  = TestAC_data.kappa_gas;
            data.LHV        = TestAC_data.LHV;
            data.inlet_pr   = TestAC_data.inlet_pr;
            data.comb_pr    = TestAC_data.comb_pr;
            data.eta_fan    = TestAC_data.eta_fan;
            data.eta_lpc    = TestAC_data.eta_lpc;
            data.eta_hpc    = TestAC_data.eta_hpc;
            data.eta_hpt    = TestAC_data.eta_hpt;
            data.eta_lpt    = TestAC_data.eta_lpt;
            data.eta_b      = TestAC_data.eta_b;
            data.eta_mech   = TestAC_data.eta_mech;
            data.nozzle_eff = TestAC_data.nozzle_eff;
            data.FAR        = TestAC_data.FAR;
            data.M21        = TestAC_data.M21;

            % -- thermodynamics outputs (filled by MDA loop) --
            data.mdot_total = NaN;
            data.TIT        = NaN;
            data.thrust     = NaN;
            data.TSFC       = NaN;
            data.fuel_flow  = NaN;
            data.A_inlet    = NaN;
            data.A_core     = NaN;

            % -- loop control --
            data.verbose    = false;   % set true to print per-discipline detail

            % -- aircraft / wing constants --
            data.S_ref   = TestAC_data.S_ref;
            data.AR      = TestAC_data.AR;
            data.e       = TestAC_data.e;
            data.sweep   = TestAC_data.sweep;
            data.t_c     = TestAC_data.t_c;
            data.taper   = TestAC_data.taper;
            data.eta_eng = TestAC_data.eta_eng;
            data.CD0_ref = TestAC_data.CD0_ref;
            data.CL_ref  = TestAC_data.CL_ref;
            data.MTOW    = TestAC_data.MTOW;
            data.OEW_ref = TestAC_data.OEW_ref;
            data.W_fuel  = TestAC_data.W_fuel;
            data.R_range      = TestAC_data.R_range;
            data.fuel_reserve = TestAC_data.fuel_reserve;
            data.n_eng        = TestAC_data.n_eng;
            data.rho_mat = TestAC_data.rho_mat;

            % -- reference engine --
            data.D_fan_ref       = TestAC_data.D_fan_ref;
            data.L_eng_ref       = TestAC_data.L_eng_ref;
            data.S_wet_ref       = TestAC_data.S_wet_ref;
            data.W_eng_ref       = TestAC_data.W_eng_ref;
            data.thrust_max_ref  = TestAC_data.thrust_max_ref;
            data.OPR_ref         = TestAC_data.OPR_ref;

            % -- constraint thresholds --
            data.TIT_max      = TestAC_data.TIT_max;
            data.M_tip_max    = TestAC_data.M_tip_max;
            data.h_engine     = TestAC_data.h_engine;
            data.clearance_min = TestAC_data.clearance_min;
            data.tip_gap_frac  = TestAC_data.tip_gap_frac;
            data.t_max        = TestAC_data.t_max;
        end

    end
end