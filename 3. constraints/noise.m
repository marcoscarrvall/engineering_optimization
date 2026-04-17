function [c_noise] = noise(state, x, atm, thermo, eng)
    % Original file used x, atm, thermo... but then we changed it
    % so that M_tip was calculated in engine sizing :)
    % ---- Relative tip Mach number --------------------------------------
    M_tip = state.M_tip;  

    state.M_tip = M_tip;

    % Normalised inequality: g(x) = (M_tip - M_tip_max) / M_tip_max <= 0
    c_noise  = (M_tip - eng.M_tip_max) / eng.M_tip_max;
    violated = c_noise > 0;

    if violated
        fprintf('  STATUS        : *** VIOLATED (fan too fast / too loud) ***\n');
    end

end