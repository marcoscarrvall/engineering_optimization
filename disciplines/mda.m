function state = mda(dv, ac, atm, eng, thermo_data, wate_data, options)
%% MDA.m  —  Multidisciplinary Analysis Converger
%
% Discipline coupling chain:
%
%   WATE  ──►  AERO  ──►  THERMO  ──┐
%    ▲                               │
%    └───────────── mdot ────────────┘
%
% Convergence variable : OEW  (Operating Empty Weight, kg)
%   OEW = W_engine_dry × N_engines  (pylon included)
%   This propagates through AERO (MTOW update) → THERMO (thrust req.) → WATE
%
% Design variable exposed : dv.BPR  (all other cycle DVs are held constant)
%
% USAGE
%   state = MDA(dv, ac, atm, eng, thermo_data, wate_data)
%   state = MDA(dv, ac, atm, eng, thermo_data, wate_data, 'conv_eps', 1e-5)
%
% INPUTS
%   dv          struct   Design variables  (BPR, PR_Fan, PR_LPC, PR_HPC, V)
%   ac          struct   Aircraft reference data
%   atm         struct   Atmospheric / ISA cruise conditions
%   eng         struct   Engine reference scalars
%   thermo_data struct   Thermodynamic constants and efficiencies
%   wate_data   struct   WATE calibration data (passed through)
%   options     name-value pairs (see arguments block below)
%
% OUTPUTS
%   state       struct   Converged state (all discipline outputs)
%
% =========================================================================

    arguments
        dv
        ac
        atm
        eng
        thermo_data
        wate_data
        options.conv_eps   (1,1) double  = 1e-4   % convergence tolerance on OEW [-]
        options.conv_niter (1,1) double  = 30      % maximum MDA iterations
        options.display    char          = 'iter'  % 'iter' | 'time' | 'iter+time' | 'none'
    end

    disp_iter = contains(options.display, 'iter');
    disp_time = contains(options.display, 'time');

    %% ── 1.  INITIALISE STATE ─────────────────────────────────────────────
    %
    % Seed values come from reference engine (CFM56-5B class).
    % Units: SI throughout.

    % ── Seed state ──────────────────────────────────────────────────────
    % mdot seed: use a rough Fsp estimate so WATE starts from a physically
    % plausible fan size rather than an arbitrary constant.
    %   Fsp_seed ≈ V_jet - V_inf.  For a high-BPR fan: V_jet ≈ V + 80 m/s.
    Fsp_seed        = 80;                            % [N/(kg/s)]  rough seed
    D_total_seed    = ac.CD0_clean * ...
                      (0.5 * atm.rho_cr * dv.V^2) * ac.S_ref * 1.6; % generous
    mdot_seed       = D_total_seed / (ac.N_engines * max(Fsp_seed, 1));

    state.mdot      = mdot_seed;       % [kg/s]  total both engines
    state.D_total   = D_total_seed;    % [N]     cruise drag seed
    state.TIT       = 1550;            % [K]
    state.TSFC      = 1.55e-5;         % [kg/N/s]
    state.CL        = ac.CL_cr;        % [-]
    state.CL_CD     = 17.0;            % [-]

    % Engine geometry / weight — seeded from reference engine
    state.W_engine  = eng.W_engine_ref;
    state.D_fan     = eng.D_fan_ref;
    state.D_nacelle = eng.D_fan_ref * 1.10;
    state.L_engine  = eng.L_engine_ref;
    state.A_fan     = pi/4 * eng.D_fan_ref^2 * (1 - 0.30^2);
    state.N_stages  = 0;
    state.MTOW      = ac.MTOW_ref;

    % Derived convergence variable
    OEW_i = ac.N_engines * state.W_engine;

    %% ── 2.  MDA LOOP ─────────────────────────────────────────────────────
    %
    % Convergence is tracked on mdot — the direct coupling variable passed
    % from THERMO into WATE.  OEW is also printed as a physically meaningful
    % secondary indicator.
    %
    %   WATE(mdot_i) → geometry/weight → AERO → drag
    %   THERMO(drag)  → mdot_{i+1}
    %   converged when |mdot_{i+1} - mdot_i| / mdot_i < conv_eps

    err      = Inf;
    iter     = 0;
    mdot_i   = state.mdot;

    t_start = tic;

    if disp_iter
        fprintf('\n════════════════════════════════════════════════════════════════\n');
        fprintf('  MDA  —  BPR = %.2f\n', dv.BPR);
        fprintf('════════════════════════════════════════════════════════════════\n');
        fprintf('  %-4s  %-12s  %-10s  %-10s  %-10s  %-10s\n', ...
                'Iter', 'mdot[kg/s]', 'Drag[N]', 'OEW[kg]', 'L/D', 'Rel.err');
        fprintf('  %s\n', repmat('-', 1, 64));
    end

    while (err > options.conv_eps) && (iter < options.conv_niter)
        iter = iter + 1;

        % ── (a) WATE: engine geometry + weight from current mdot ──────────
        state = wate(dv, atm, eng, state, false);

        % ── (b) AERO: drag + L/D from updated engine weight & nacelle size ─
        state = aero(state, ac, atm, dv, eng, false);

        % ── (c) THERMO: find mdot so that F_net = D_total/2 ───────────────
        %   state.mdot is updated → feeds back into WATE on next iteration
        state = thermo(state, dv, atm, thermo_data, false);

        % ── Convergence check on mdot ─────────────────────────────────────
        mdot_j = state.mdot;
        err    = abs(mdot_j - mdot_i) / max(abs(mdot_i), 1);
        OEW_j  = ac.N_engines * state.W_engine;

        if disp_iter
            fprintf('  %-4d  %-12.3f  %-10.1f  %-10.1f  %-10.4f  %-10.3e\n', ...
                    iter, mdot_j, state.D_total, OEW_j, state.CL_CD, err);
        end

        mdot_i = mdot_j;
        OEW_i  = OEW_j;
    end

    t_conv = toc(t_start);

    %% ── 3.  CONVERGENCE REPORT ───────────────────────────────────────────
    if iter >= options.conv_niter && err > options.conv_eps
        warning('MDA:NotConverged', ...
            'MDA did not converge after %d iterations (err = %.3e).', iter, err);
    end

    if disp_iter || disp_time
        fprintf('\n');
    end
    if disp_iter
        status = 'CONVERGED';
        if iter >= options.conv_niter && err > options.conv_eps
            status = 'NOT CONVERGED';
        end
        fprintf('  Status  : %s\n', status);
        fprintf('  Iters   : %d / %d\n', iter, options.conv_niter);
        fprintf('  Final err: %.3e  (tol = %.3e)\n', err, options.conv_eps);
    end
    if disp_time
        fprintf('  Wall time: %.2f s\n', t_conv);
    end
    if disp_iter || disp_time
        fprintf('════════════════════════════════════════════════════════\n\n');
    end

    %% ── 4.  FINAL VERBOSE DISCIPLINE PASS ───────────────────────────────
    state = wate(dv, atm, eng, state, true);
    state = aero(state, ac, atm, dv, eng, true);
    state = thermo(state, dv, atm, thermo_data, true);

    %% ── 5.  STORE LOOP METADATA ──────────────────────────────────────────
    state.mda.BPR        = dv.BPR;
    state.mda.iter       = iter;
    state.mda.err        = err;
    state.mda.converged  = (err <= options.conv_eps);
    state.mda.wall_time  = t_conv;
    state.mda.OEW        = OEW_i;            % converged OEW

end