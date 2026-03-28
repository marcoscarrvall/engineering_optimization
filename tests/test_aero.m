%% test_aero.m
clear; clc;
data = TestAC_data;

x.BPR = 10.0;
x.V   = sqrt(data.atm.gamma * data.atm.R * data.atm.T_cruise) * 0.82;  

%% ---- 2. Pass-through structs -------------------------------------------
atm = data.atm;
eng = data.eng;
ac  = data.ac;

%% ---- 3. Build a realistic state ----------------------------------------
% Plug in values close to a real narrow-body (A320-class) at cruise
state           = data.state;
state.D_fan     = 1.80;        % [m]    slightly larger than ref (1.73 m)
state.L_eng     = 3.60;        % [m]    slightly longer than ref (3.50 m)
state.S         = 122.6;       % [m^2]  wing area = reference

%% ---- 4. Call aero ------------------------------------------------------
state = aero(state, x, atm, eng, ac, true);
