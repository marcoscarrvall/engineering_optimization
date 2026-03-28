clear; clc;
data = TestAC_data;

x.BPR = 6.0;
x.V   = sqrt(data.atm.gamma * data.atm.R * data.atm.T_cruise) * 0.78;

x_consts.PR_fan = 1.7;
x_consts.PR_LPC = 2.6;
x_consts.PR_HPC = 6.1;

options.tol      = 1e-4;
options.max_iter = 50;
options.verbose  = true;

state = mda(x, x_consts, data, options);

%% ---- Summary -----------------------------------------------------------
fprintf('=== CONVERGED STATE ===\n');
fprintf('  MTOW        = %8.1f kg\n',    state.MTOW);
fprintf('  D_cruise    = %8.1f N\n',     state.D_cruise);
fprintf('  mdot        = %8.2f kg/s\n',  state.mdot);
fprintf('  TIT         = %8.1f K\n',     state.TIT);
fprintf('  TSFC        = %.4e kg/N/s\n', state.TSFC);
fprintf('  D_fan       = %8.4f m\n',     state.D_fan);
fprintf('  L_eng       = %8.4f m\n',     state.L_eng);
fprintf('  N_stages    = %8d\n',         state.N_stages);
fprintf('  W_engine    = %8.1f kg\n',    state.W_engine);
fprintf('  W_wing      = %8.1f kg\n',    state.W_wing);
fprintf('  L/D         = %8.4f\n',       state.CL_CD);
fprintf('  S           = %8.2f m^2\n',   state.S);