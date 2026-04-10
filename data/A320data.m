classdef A320data
    % A320DATA Configuration constants for Airbus A320 Engine/Aircraft Analysis
    
    properties (Constant)
        name = "A320"

        atm = struct( ...
            "rho_cruise",   0.3639,     ... % [kg/m^3]  ISA density
            "T_cruise",     218.8,      ... % [K]       ISA static temperature
            "P_cruise",     23840.0,    ... % [Pa]      ISA static pressure
            "R",            287.05      ... % [J/kg/K]  Specific gas constant, air
        )

        thermo_data = struct( ...
            "Cp_air",       1005,       ... % [J/kgK]   Cp cold side
            "Cp_gas",       1150,       ... % [J/kgK]   Cp hot side
            "gamma_c",      1.4,        ... % [-]       Gamma cold
            "gamma_h",      1.33,       ... % [-]       Gamma hot
            "LHV",          43e6,       ... % [J/kg]    Lower Heating Value (Jet A-1)
            "eta_cc",       0.99,       ... % [-]       Combustion efficiency
            "dP_cc_frac",   0.04,       ... % [-]       Burner pressure drop fraction
            "eta_LPC",      0.89,       ... % [-]       LPC isentropic efficiency
            "eta_HPC",      0.90,       ... % [-]       HPC isentropic efficiency
            "eta_HPT",      0.92,       ... % [-]       HPT isentropic efficiency
            "eta_LPT",      0.92,       ... % [-]       LPT isentropic efficiency
            "eta_mech",     0.98,       ... % [-]       Mechanical efficiency
            "FAR",          0.025       ... % [-]       Fuel-Air Ratio
        )

        ac = struct( ...
            "N_eng",            2,          ... % [-]       Number of engines
            "W_fuel",           19000,      ... % [kg]      Max usable fuel
            "OEW_ref",          40000,      ... % [kg]      Operating Empty Weight
            "W_pay",            17600,      ... % [kg]      Design payload
            "S_ref",            122.6,      ... % [m^2]     Wing reference area
            "AR",               9.39,       ... % [-]       Aspect Ratio
            "e",                0.85,       ... % [-]       Oswald efficiency
            "CL_cruise_ref",    0.649,      ... % [-]       Ref Cruise Lift Coeff
            "CD_cruise_ref",    0.038,      ... % [-]       Ref Cruise Drag Coeff
            "CD_parasitic_ref", 0.00171,    ... % [-]       Ref Nacelle parasitic drag
            "range",            9192.1e3    ... % [m]       Design range
        )

        eng = struct( ...
            "W_eng_ref",    3000,       ... % [kg]      Reference engine dry mass
            "D_fan_ref",    1.933,      ... % [m]       Reference fan diameter
            "L_eng_ref",    3.728,      ... % [m]       Reference engine length
            "rho_mat",      4430,       ... % [kg/m^3]  Density of Titanium
            "eta_fan",      0.91        ... % [-]       Fan isentropic efficiency
        )

        constraints = struct( ...
            "h_engine",         1.4, ... % [m] Engine center line position above ground
            "min_clearance",    0.35, ... % [m] Minimum ground clearance
            "m_tip_max",        1.55,   ... % [-] Maximum fan tip Mach number (to avoid shock waves)
            "TIT_max",        1800   ... % [K] Maximum Turbine Inlet Temperature
        )

        state = struct( ...
            "eta_fan",  0.91,   ... % [-]       Initial fan efficiency
            "D_cruise", 44103.5   ... % [N]       Initial cruise thrust required
        )
    end
end