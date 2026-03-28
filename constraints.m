function [g, h] = constraints(x)
    % Perform MDA again (or use a global variable/cache to save time)
    y = solve_mda(x);
    
    % 3 Constraint Functions (g(x) <= 0 format)
    g = zeros(3, 1);
    g(1) = y(1) - 50;   % Example: Limit on a discipline output
    g(2) = x(2) + y(2); % Example: Structural constraint
    g(3) = 1 - y(4);    % Example: Minimum performance requirement
    
    h = []; % No equality constraints required for MDF structure
end