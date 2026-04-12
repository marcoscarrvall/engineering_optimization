function objective = optim(x, lb, ub, data, mda_options)

    x = denormalize_vars(x, lb, ub);
    global optHistory

    state_converged = mda(x, data, mda_options);
    
    objective = real(double((data.ac.range - state_converged.range) / data.ac.range));
    
    optHistory.fval(end+1) = objective;
    optHistory.x(end+1,:) = x;
    optHistory.iter(end+1) = length(optHistory.fval);


    function x_real = denormalize_vars(x_norm, lb, ub)
        % Scales [0, 1] values back to physical units (for the MDA)
        x_real = x_norm .* (ub - lb) + lb;
    end
end 