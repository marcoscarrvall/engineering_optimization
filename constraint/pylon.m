function [c_pylon] = pylon(Total_Weight, thrust_max)
    % Inputs: 
    % Total_Weight: From Discipline 3 [kg]
    % thrust_max: Max Takeoff Thrust [N]

    % Load factor for structural safety (standard maneuver limit + safety margin)
    n_z = 4.5; 
    
    % K factor for a conventional under-wing pylon (A320 class)
    K_pyl = 0.065; 

    % Pylon weight scaling law
    % The weight increases non-linearly with the mass it must carry
    W_pylon = K_pyl * (Total_Weight * n_z)^0.78 * (thrust_max^0.12);
    
    % Direct Peer Note: If you increase OPR, engine weight goes up, 
    % which cascades into a heavier pylon, increasing total aircraft OEi.
    %% 4. Pylon Constraint Check (A320 Integration)
    % Max allowable weight for the A320 CFM56/LEAP class pylon
    c_pylon = 2500 - W_pylon;
    
    % Direct Peer Note: If Margin < 0, your OPR/BPR optimization 
    % has made the engine physically too large for the wing structure.
end