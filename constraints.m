function [g, h] = constraints(x, TestAC_data)

    state_converged = mda(x);
    g = zeros(3, 1);
    g(1) = clearance(state_converged, TestAC_data.eng, TestAC_data.con);
    g(2) = noise(state_converged, TestAC_data.dv, TestAC_data.atm, TestAC_data.thermo, TestAC_data.con);
    g(3) = tit(state_converged, TestAC_data.con);
    
    h = []; 
end