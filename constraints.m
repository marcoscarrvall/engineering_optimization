function [g, h] = constraints(x, TestAC_data, options_mda)

    state_converged = mda(x, TestAC_data, options_mda);
    g = zeros(3, 1);
    g(1) = clearance(state_converged, TestAC_data.eng);
    g(2) = noise(state_converged, x, TestAC_data.atm, TestAC_data.thermo, TestAC_data.eng);
    g(3) = tit(state_converged, TestAC_data.eng);
    
    h = []; 
end