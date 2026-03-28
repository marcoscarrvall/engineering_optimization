function range = optim(x, x_consts, TestAC_data, options_mda)

    state_converged = mda(x, x_consts, TestAC_data, options_mda);
    
    state_converged = breguet(state_converged, x, TestAC_data.atm, TestAC_data.ac, TestAC_data.mission);
    fprintf('  Range = %.1f km\n', state_converged.range);
    range = 1-(state_converged.range/2346747.3);
end 