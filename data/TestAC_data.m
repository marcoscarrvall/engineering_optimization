classdef TestAC_data

    properties (Constant)

        name = "A320"

        atm = struct( ...
            "rho_cruise",   0.3639,    ... % [kg/m^3]  ISA density at 35,000 ft
            "T_cruise",     218.8,     ... % [K]       ISA static temperature at 35,000 ft
            "P_cruise",     23842.0,   ... % [Pa]      ISA static pressure at 35,000 ft
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
            "LHV",          43e6,      ... % [J/kg]    fuel lower heating value (Jet-A)
            "FAR",          0.02755,   ... % [-]       fuel-to-air ratio (~900 K temp rise, CFM56)
            ... % --- Combustor ---
            "eta_cc",       0.995,     ... % [-]  combustor efficiency
            "dP_cc_frac",   0.04,      ... % [-]  combustor total-pressure loss
            "inlet_pr",     0.995,     ... % [-]  inlet total-pressure recovery
            ... % --- Component isentropic efficiencies (CFM56-5B level) ---
            "eta_fan",      0.92,      ... % [-]  fan
            "eta_LPC",      0.91,      ... % [-]  low-pressure compressor (3 stages)
            "eta_HPC",      0.90,      ... % [-]  high-pressure compressor (9 stages)
            "eta_HPT",      0.91,      ... % [-]  high-pressure turbine
            "eta_LPT",      0.92,      ... % [-]  low-pressure turbine
            "eta_mech",     0.995,     ... % [-]  mechanical (shaft)
            "nozzle_eff",   0.99,      ... % [-]  nozzle isentropic efficiency
            ... % --- Mach at LPC face ---
            "M21",          0.55       ... % [-]  Mach number at station 21 (LPC face)
        )

        eng = struct( ...
            ... % CFM56-5B reference engine
            "D_fan_ref",        1.735,     ... % [m]    CFM56-5B fan diameter
            "L_eng_ref",        2.60,      ... % [m]    CFM56-5B overall length
            "S_wet_ref",        14.2,      ... % [m^2]  nacelle wetted area (≈ pi·D·L)
            "W_eng_ref",        2630,      ... % [kg]   CFM56-5B installed + pylon mass
            "thrust_max_ref",   120000,    ... % [N]    CFM56-5B max take-off thrust (27,000 lbf)
            "OPR_ref",          27.0,      ... % [-]    CFM56-5B overall pressure ratio
            "rho_mat",          4430,      ... % [kg/m^3]  titanium alloy
            "TIT_max",          1800.0,    ... % [K]    max allowable TIT
            "M_tip_max",        1.05,      ... % [-]    max fan tip relative Mach number
            "h_engine",         2.10,      ... % [m]    engine centreline height above ground
            "clearance_min",    0.30,      ... % [m]    minimum fan tip-to-ground clearance
            "tip_gap_frac",     0.005      ... % [-]    tip-to-nacelle gap fraction
        )

        wate = struct( ...
            "hub_to_tip",       0.30,      ... % [-]    fan hub-to-tip radius ratio
            "M_axial",          0.55,      ... % [-]    axial Mach number at fan face
            "p_stage",          0.20,      ... % [m]    axial length per compressor stage
            "L_combustor",      0.35,      ... % [m]    combustor axial length
            "W_misc",           150        ... % [kg]   miscellaneous systems mass
        )

        ac = struct( ...
            ... % --- Wing geometry (A320ceo) ---
            "S_ref",        135,     ... % [m^2]  wing reference area
            "AR",           9.39,      ... % [-]    wing aspect ratio
            "e",            0.85,      ... % [-]    Oswald span efficiency
            "sweep",        25.0,      ... % [deg]  wing quarter-chord sweep
            "t_c",          0.115,     ... % [-]    wing thickness-to-chord ratio
            "taper",        0.24,      ... % [-]    wing taper ratio
            "eta_eng",      0.34,      ... % [-]    engine spanwise station (2y/b)
            ... % --- Aerodynamics ---
            "CL_cruise",    0.50,      ... % [-]    cruise CL  (M0.78, 35 000 ft, mid-cruise weight)
            "CD_cruise",    0.0285,    ... % [-]    cruise CD  (L/D ≈ 17.5, consistent with CL)
            "CD_parasitic_ref", 0.0160,... % [-]    zero-lift CD baseline
            ... % --- Weights (A320ceo, MTOW variant) ---
            "MTOW",         77000,     ... % [kg]   maximum take-off weight
            "OEW",          42600,     ... % [kg]   operating empty weight
            "W_fuel",       18000,     ... % [kg]   max usable fuel (A320ceo tanks)
            "W_pay",        16500,     ... % [kg]   design payload (150 pax × 110 kg)
            "W_wing",       9500,      ... % [kg]   wing structural mass
            ... % --- Propulsion ---
            "N_eng",        2          ... % [-]    number of engines
        )

        mission = struct( ...
            ... % --- Segment fuel fractions ---
            "ff_takeoff",   0.970,     ... % [-]  take-off
            "ff_climb",     0.985,     ... % [-]  climb to cruise altitude
            "ff_descent",   0.992,     ... % [-]  descent and landing
            ... % --- Reserve ---
            "fuel_reserve", 0.05,      ... % [-]  5 % reserve fuel
            ... % --- Design mission ---
            "R_range",      3300e3,    ... % [m]   A320 design range (3 300 km / 1 780 NM)
            ... % --- Endurance ---
            "t_max",        8.0*3600   ... % [s]   max flight time (~8 h for A320)
        )

        state = struct( ...
            ... % --- Mass flows ---
            "m_dot_total",  NaN,   ... % [kg/s]  total engine mass-flow (fan face)
            ... % --- Weights (sized) ---
            "MTOW",         NaN,   ... % [kg]   max take-off weight (converged)
            "W_eng",        NaN,   ... % [kg]   installed engine weight, per engine
            "W_wing",       NaN,   ... % [kg]   wing structural weight
            "S",            NaN,   ... % [m^2]  wing area
            ... % --- Engine geometry (sized) ---
            "D_fan",        NaN,   ... % [m]    fan tip diameter
            "D_nacelle",    NaN,   ... % [m]    nacelle outer diameter
            "L_eng",        NaN,   ... % [m]    overall engine length
            "A_fan",        NaN,   ... % [m^2]  fan annulus area
            ... % --- Compressor stages ---
            "N_stages",     NaN,   ... % [-]    total compressor stage count
            ... % --- Cruise aerodynamics ---
            "D_cruise",     NaN,   ... % [N]    cruise drag
            "CL",           NaN,   ... % [-]    cruise lift coefficient
            "CD",           NaN,   ... % [-]    cruise drag coefficient
            "CL_CD",        NaN,   ... % [-]    lift-to-drag ratio
            ... % --- Performance ---
            "TIT",          NaN,   ... % [K]    turbine inlet temperature
            "TSFC",         NaN,   ... % [kg/N/s]  thrust-specific fuel consumption
            "range",        NaN,   ... % [m]    achievable range (Breguet)
            "clearance",    NaN    ... % [m]    fan tip-to-ground clearance
             )

    end
end