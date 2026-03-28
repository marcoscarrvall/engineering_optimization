function [c_cl] = clearance(state, eng)

    clearance_val = eng.h_engine - state.D_nacelle;

    state.clearance = clearance_val;  

    c_cl     = (eng.clearance_min - clearance_val) / eng.clearance_min;
end