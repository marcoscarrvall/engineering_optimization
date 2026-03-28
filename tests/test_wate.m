clear; clc;
data = TestAC_data;

x.BPR = 5.0;
x.V   = 230;          

x_consts.PR_fan  = 1.6;
x_consts.PR_LPC  = 3.0;
x_consts.PR_HPC  = 5.6;

thermo_data = data.thermo;
atm = data.atm;
ac = data.ac;
eng = data.eng;
wate_consts = data.wate;

state        = data.state;
state.D_total = 25000;


state = thermo(state, x, x_consts, atm, thermo_data, ac, false);

state = wate(state, x, x_consts, atm, eng, wate_consts, ac, true);

%% ---- 10. Summary -------------------------------------------------------
fprintf('\n=== SUMMARY ===\n');
fprintf('  BPR             = %.1f\n',   x.BPR);
fprintf('  mdot total      = %.2f kg/s\n', state.mdot);
fprintf('  D_fan           = %.3f m\n', state.D_fan);
fprintf('  L_engine        = %.3f m\n', state.L_engine);
fprintf('  N_stages        = %d\n',     state.N_stages);
fprintf('  W_engine        = %.1f kg\n', state.W_engine);
fprintf('  MTOW            = %.1f kg\n', state.MTOW);
fprintf('  W_wing          = %.1f kg\n', state.W_wing);