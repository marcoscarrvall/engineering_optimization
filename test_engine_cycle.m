% test_engine_cycle.m
% Tests engine_cycle.m and breguet_range.m across a range of BPR and V0.
% Fuel weight is fixed; output is the achievable range.

clear; clc;

d = data();

n_engines = 2;
thrust_N  = d.W_start / d.CL_CD / n_engines;

fprintf('=================================================================\n');
fprintf('  engine_cycle.m + breguet_range.m — test cases\n');
fprintf('  Fuel weight: %.0f kN  |  L/D: %.1f  |  W_start: %.0f kN\n', ...
        d.W_fuel/1e3, d.CL_CD, d.W_start/1e3);
fprintf('  Thrust per engine (cruise): %.1f kN\n', thrust_N/1e3);
fprintf('=================================================================\n');
fprintf('  %-5s  %-8s  |  %-8s  %-8s  %-8s  |  %-10s\n', ...
        'BPR', 'V0(m/s)', 'comp_pr', 'TIT(K)', 'TSFC', 'Range(km)');
fprintf('  %s\n', repmat('-', 1, 63));

cases = [
    %  BPR    V0
     6.0,   230;
    10.0,   230;
    12.5,   230;
    15.0,   230;
    12.5,   180;
    12.5,   230;
    12.5,   270;
];

for i = 1:size(cases, 1)
    BPR = cases(i, 1);
    V0  = cases(i, 2);
    design_vec = [BPR, V0];

    try
        [comp_pr, TIT, TSFC] = engine_cycle(d, design_vec, thrust_N);
        range_km = breguet_range(d, design_vec, d.W_start, thrust_N) / 1e3;

        fprintf('  %-5.1f  %-8.1f  |  %-8.4f  %-8.1f  %-8.4f  |  %-10.1f\n', ...
                BPR, V0, comp_pr, TIT, TSFC, range_km);
    catch ME
        fprintf('  %-5.1f  %-8.1f  |  ERROR: %s\n', BPR, V0, ME.message);
    end

    if i == 4
        fprintf('  %s\n', repmat('.', 1, 63));
    end
end

fprintf('=================================================================\n');