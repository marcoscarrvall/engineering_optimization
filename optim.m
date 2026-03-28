function f = optim(x)
    % Perform MDA to ensure multidisciplinary feasibility
    y = run_mda(x);
    
    % Define your objective here (e.g., minimize weight or maximize efficiency)
    % f = ... your formula using x and y ...
    f = x(1)^2 + y(3); % Placeholder
end