function [g, h] = constraints(x, lb, ub, data, mda_options)
    global optHistory

    
    x = denormalize_vars(x, lb, ub);

    state = mda(x, data, mda_options);

    g = zeros(3, 1);
    clearance = data.constraints.h_engine - state.D_fan/2;
    g(1) = real(double((data.constraints.min_clearance - clearance) / data.constraints.min_clearance));
    g(2) = real(double((-data.constraints.TIT_max + state.TIT) / data.constraints.TIT_max));
    g(3) = real(double((-data.constraints.m_tip_max + state.M_tip) / data.constraints.m_tip_max));


    if ~isempty(optHistory.fval)
        curr_idx = length(optHistory.fval);
        % Store the constraint values corresponding to the latest objective evaluation
        optHistory.constr(curr_idx, :) = g';
    end

    h = []; 
    function x_real = denormalize_vars(x_norm, lb, ub)
        % Scales [0, 1] values back to physical units (for the MDA)
        x_real = x_norm .* (ub - lb) + lb;
    end
end