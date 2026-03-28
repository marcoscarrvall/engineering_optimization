function range = optim(x, TestAC_data)
    state_converged = run_mda(x);
    
    state_converged = breguet(state_converged, TestAC_data.dv, TestAC_data.ac, TestAC_data.atm, TestAC_data.mis, TestAC_data.print_flag);

    range = state_converged.R;
end