classdef TestAC_data

    properties (Constant)

        name = "testAC"

        atm = struct( ...
            "rho_cruise",   0.3119,    ... % [kg/m^3]  cruise air density
            "T_cruise",     216.5,     ... % [K]       ISA tropopause static temperature
            "P_cruise",     22632.0,   ... % [Pa]      ISA tropopause static pressure
            "Mach_cr",      0.78,      ... % [-]       cruise Mach number
            "V0",           230.053327,... % [m/s]     freestream velocity
            "gamma",        1.4,       ... % [-]       cp/cv, air
            "R",            287.0,     ... % [J/kg/K]  specific gas constant, air
            "g",            9.81       ... % [m/s^2]   gravitational acceleration
        )

        thermo = struct( ...
            ... % --- Gas properties ---
            "Cp_air",       1005.0,    ... % [J/kg/K]  specific heat, cold stream (air)
            "Cp_gas",       1150.0,    ... % [J/kg/K]  specific heat, hot stream (combustion gas)
            "gamma_c",      1.4,       ... % [-]       cp/cv, cold stream
            "gamma_h",      1.33,      ... % [-]       cp/cv, hot stream
            "LHV",          43e6,      ... % [J/kg]    fuel lower heating value
            "FAR",          0.023927, ...% [-]     fuel-to-air ratio (pre-computed, delta_T = 900 K)
            ... % --- Combustor ---
            "eta_cc",       0.995,     ... % [-]  combustor efficiency          (= eta_b)
            "dP_cc_frac",   0.04,      ... % [-]  combustor total-pressure loss (= 1 - comb_pr)
            "inlet_pr",     0.99,      ... % [-]  inlet total-pressure recovery
            ... % --- Component isentropic efficiencies ---
            "eta_fan",      0.92,      ... % [-]  fan
            "eta_LPC",      0.92,      ... % [-]  low-pressure compressor
            "eta_HPC",      0.92,      ... % [-]  high-pressure compressor
            "eta_HPT",      0.91,      ... % [-]  high-pressure turbine
            "eta_LPT",      0.91,      ... % [-]  low-pressure turbine
            "eta_mech",     0.995,     ... % [-]  mechanical (shaft)
            "nozzle_eff",   0.99,      ... % [-]  nozzle isentropic efficiency
            ... % --- Mach at LPC face (area sizing) ---
            "M21",          0.8        ... % [-]  Mach number at station 21 (LPC face)
        )

        eng = struct( ...
            "D_fan_ref",        1.73,      ... % [m]    fan tip diameter
            "L_eng_ref",        3.50,      ... % [m]    overall engine length
            "S_wet_ref",        19.0,      ... % [m^2]  nacelle wetted area (≈ pi·D·L)
            "W_eng_ref",        2400,      ... % [kg]   installed engine + pylon, per engine
            "thrust_max_ref",   120000,    ... % [N]    max take-off thrust, per engine
            "OPR_ref",          27.0,      ... % [-]    overall pressure ratio, reference engine
            "rho_mat",          4430,      ... % [kg/m^3]  titanium alloy (fan/compressor discs)
            "TIT_max",          1800.0,    ... % [K]    max allowable turbine inlet temperature
            "M_tip_max",        1.05,      ... % [-]    max fan tip relative Mach number
            "h_engine",         2.10,      ... % [m]    engine centreline height above ground
            "clearance_min",    0.30,      ... % [m]    minimum fan tip-to-ground clearance
            "tip_gap_frac",     0.005      ... % [-]    min tip-to-nacelle gap as fraction of D_fan
        )

        wate = struct( ...
            "hub_to_tip",       0.30,      ... % [-]    fan hub-to-tip radius ratio
            "M_axial",          0.50,      ... % [-]    axial Mach number at compressor face
            "p_stage",          1.35,      ... % [-]    pressure ratio per fan stage
            "L_combustor",      0.35,      ... % [m]    combustor length (annular)
            "W_misc",           150        ... % [kg]   miscellaneous engine systems mass
        )

        ac = struct( ...
            ... % --- Wing geometry ---
            "S_ref",        122.6,     ... % [m^2]  wing reference area
            "AR",           9.5,       ... % [-]    wing aspect ratio
            "e",            0.85,      ... % [-]    Oswald span efficiency
            "sweep",        25.0,      ... % [deg]  wing quarter-chord sweep
            "t_c",          0.12,      ... % [-]    wing thickness-to-chord ratio
            "taper",        0.25,      ... % [-]    wing taper ratio
            "eta_eng",      0.35,      ... % [-]    engine spanwise station (2y/b)
            ... % --- Aerodynamics ---
            "CL_cruise",    0.50,      ... % [-]    cruise lift coefficient
            "CD0_ref",      0.0200,    ... % [-]    baseline zero-lift drag coefficient
            ... % --- Weights ---
            "MTOW",         77000,     ... % [kg]   max take-off weight
            "OEW",          42000,     ... % [kg]   operating empty weight (baseline engine)
            "W_fuel",       20000,     ... % [kg]   fuel load at start of cruise
            "W_pay",        18000,     ... % [kg]   design payload  ← set representative value
            "W_wing",       8500,      ... % [kg]   wing structural mass ← set representative value
            ... % --- Propulsion layout ---
            "N_eng",        2          ... % [-]    number of engines
        )

        mission = struct( ...
            ... % --- Segment fuel fractions (W_end / W_start per segment) ---
            "ff_takeoff",   0.970,     ... % [-]  take-off   (engine start → brake release)
            "ff_climb",     0.985,     ... % [-]  climb      (brake release → TOC)
            "ff_descent",   0.990,     ... % [-]  descent    (TOD → landing)
            ... % --- Reserve ---
            "fuel_reserve", 0.05,      ... % [-]  reserve fuel fraction (5 % of block fuel)
            ... % --- Range ---
            "R_range",      5500e3,    ... % [m]   design range (5 500 km)
            ... % --- Endurance limit ---
            "t_max",        14.0*3600  ... % [s]   max flight time (14 h crew duty limit)
        )

        state = struct( ...
            ... % --- Mass flows ---
            "m_dot_total",  NaN,   ... % [kg/s]  total engine mass-flow (fan face)
            ... % --- Weights (sized) ---
            "MTOW",         NaN,   ... % [kg]   max take-off weight (converged)
            "W_eng",        NaN,   ... % [kg]   installed engine weight, per engine
            "W_wing",       NaN,   ... % [kg]   wing structural weight
            ... % --- Engine geometry (sized) ---
            "D_fan",        NaN,   ... % [m]    fan tip diameter
            "D_nacelle",    NaN,   ... % [m]    nacelle outer diameter
            "L_engine",     NaN,   ... % [m]    overall engine length
            "A_fan",        NaN,   ... % [m^2]  fan annulus area
            ... % --- Compressor stages ---
            "N_stages",     NaN,   ... % [-]    total compressor stage count
            ... % --- Cruise aerodynamics ---
            "D_cruise",     NaN,   ... % [N]    cruise drag
            "CL",           NaN,   ... % [-]    cruise lift coefficient
            "CD",           NaN,   ... % [-]    cruise drag coefficient
            "LoD",          NaN,   ... % [-]    lift-to-drag ratio (L/D)
            ... % --- Performance ---
            "TIT",          NaN,   ... % [K]    turbine inlet temperature
            "TSFC",         NaN,   ... % [kg/N/s]  thrust-specific fuel consumption
            "range",        NaN    ... % [m]    achievable range (Breguet)
        )

    end
end