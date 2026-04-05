function range = optim(x, TestAC_data, options_mda)

    state_converged = mda(x, TestAC_data, options_mda);
    
    state_converged = breguet(state_converged, x, TestAC_data.atm, TestAC_data.ac, TestAC_data.mission);
    range = 1-(state_converged.range/10153468.473336);
end 