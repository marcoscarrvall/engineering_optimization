
classdef TestAC

    properties (Constant)
        name = "testAC"
        % wing material
        material = struct("E_modulus", 7.1e10, ...  % [Pa]
                          "density", 2800, ...      % [Kg/m3]
                          "tens_stress", 4.8e8, ... % [Pa]
                          "comp_stress", 4.6e8)     % [Pa]
        
        % airfoil file
        af_fpath = "data/BACXXX.dat"
        nCST    = 6
        % wing
        AR = 6.61   % [-]
        b  = 29   % [m]
        sweep_LE = 38.3 % [deg]
        taper  = 0.255 % [-]
        b0 = 4.06 % [m] 14% of the span
        sweep_kink = 0.01 % [deg]
        dihedral = -3 % [deg]
        twist_k = -1.2 % [deg]
        twist_t = -3 % [deg]
        incidence = 1 % [deg]

        spar_LE = 0.2 % [x/c]
        spar_TE = 0.8 % [x/c]
        panel_eff = 0.96 % [-]
        rib_pitch = 0.5 % [m]

        tank_eta0 = 0.0 % [2y/b]
        tank_eta1 = 0.7 % [2y/b]

        f_tank = 0.93 % [-]
        rho_fuel = 0.81715e3 % [Kg/m3]

        CL_CD = NaN % [NAN-> compute with range]
        
        engine = struct("weight", 1600, ... % [Kg] Weight 
                        "mounted_on", "wing", ... % option for engine placement ("wing", "fuselage")
                        "eta", 0.25, ... % [y2/b] wing position (NaN for non-wing)
                        "C_T", 1.8639e-4 ... % [-] thrust coefficient
                        ) 
        
        M_cr = 0.77 % [-] 
        h_cr = 9000 % [m]
        V_MO = 263.89 % [m/s]
        n_max = 2.5 % [-]

        MTOM = 47000 % [Kg]
        W_fuel = 10000 % [Kg]
        fuel_frac = 0.938 % [-]
        R = 3200000 % [m]
    end
end