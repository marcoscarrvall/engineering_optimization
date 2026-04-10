function objective = optim(x, data, mda_options)
    state_converged = mda(x, data, mda_options);
    fprintf("\nConverged range: %.2f km", state_converged.range/1000);
    fprintf("\nOriginal range: %.2f km", data.ac.range/1000);
    objective = (data.ac.range-state_converged.range) / data.ac.range;
    fprintf("\nObjective value: %.4f\n", objective);
end 