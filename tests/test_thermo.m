clear; clc;
data = TestAC_data;

%% ---- 1. User inputs ---------------------------------------
BPR  = 12.5;          % bypass ratio to test      [-]
D   = 30000; % cruise drag (total) to test [N]
V    = 200;     % flight velocity to test   [m/s]


%% ---- 3. Build x_consts (fixed pressure ratios) -------------------------
x_consts.PR_fan  = 1.35;    % typical turbofan values
x_consts.PR_LPC  = 4.5;
x_consts.PR_HPC  = 5.5;

x.BPR = BPR;
x.V   = V;
state.D_total = D;



state = thermo(state, x, x_consts, data.atm, data.thermo, data.ac, true);

fprintf('%-6s  %-10s  %-8s  |  %-10s  %-12s  %-10s\n', ...
        'BPR', 'D [N]', 'V [m/s]', 'mdot [kg/s]', 'TSFC [mg/N/s]', 'TIT [K]');
fprintf('%s\n', repmat('-', 1, 72));
fprintf('%-6.1f  %-10.0f  %-8.1f  |  %-10.2f  %-12.4f  %-10.1f\n', ...
        BPR, D, V, state.mdot, state.TSFC * 1e6, state.TIT);

