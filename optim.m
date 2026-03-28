function range = optim(x, x_consts, TestAC_data, options_mda)

    state_converged = mda(x, x_consts, TestAC_data, options_mda);
    
    state_converged = breguet(state_converged, x, TestAC_data.atm, TestAC_data.ac, TestAC_data.mission);

    range = state_converged.range;
end