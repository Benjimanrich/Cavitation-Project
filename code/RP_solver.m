%% Part B: Experimental Comparison of an Acoustically Driven Bubble
% Fast, native-speed solver for both GNU Octave (lsode) and MATLAB (ode15s)
clear; clc; close all;

% 1. Load Experimental Dataset
csv_path = 'experimental_data_fig4.csv';
fid = fopen(csv_path, 'r');
if fid == -1
    error('Cannot find %s!', csv_path);
end
header_line = fgetl(fid);
fclose(fid);

if ~isempty(strfind(header_line, 'Label'))
    raw_exp = dlmread(csv_path, ',', 1, 2);
else
    raw_exp = dlmread(csv_path, ',', 1, 0);
end
t_exp = raw_exp(:, 1);    % Experimental shifted time [us]
R_exp = raw_exp(:, 2);    % Experimental bubble radius [um]
fprintf('Loaded %d experimental data points.\n', length(t_exp));

% 2. Governing Parameters in Base SI Units
R0    = 1.72e-6;          % Equilibrium bubble radius [m]
f     = 87.6e3;           % Acoustic driving frequency [Hz]
PA    = 1.72e5;           % Acoustic pressure amplitude [Pa]
gamma = 1.65;             % Polytropic exponent [-]
p0    = 1.013e5;          % Ambient static pressure [Pa]
rho   = 1000.0;           % Liquid density [kg/m^3]
sigma = 0.072;            % Surface tension [N/m]
mu    = 1.0e-3;           % Dynamic viscosity [Pa*s]
pv    = 2.3e3;            % Vapor pressure [Pa]
c     = 1500.0;           % Speed of sound in water [m/s]

pg0 = p0 + 2*sigma/R0 - pv;   % Equilibrium gas pressure [Pa]
T   = 1 / f;                  % Forcing period [s]

% 3. Integrate Stiff Rayleigh-Plesset ODE (Native Fortran lsode in Octave / ode15s in MATLAB)
n_cycles = 3;
t = linspace(0, n_cycles * T, 60000)';   % Dense uniform grid (60,000 points)

if exist('OCTAVE_VERSION', 'builtin')
    fflush(stdout);
    % Note: lsode takes f(y, t) with state vector y FIRST, time t SECOND
    rp_lsode = @(y, t_val) [ ...
        y(2); ...
        ((pv + pg0 * (R0 / y(1))^(3 * gamma) - p0 - PA * sin(2 * pi * f * t_val) ...
          - 2 * sigma / y(1) - 4 * mu * y(2) / y(1)) / rho - 1.5 * y(2)^2) / y(1) ...
    ];
    lsode_options('relative tolerance', 1e-9);
    lsode_options('absolute tolerance', 1e-12);
    y = lsode(rp_lsode, [R0; 0.0], t);
else
    rp_matlab = @(t_val, y) [ ...
        y(2); ...
        ((pv + pg0 * (R0 / y(1))^(3 * gamma) - p0 - PA * sin(2 * pi * f * t_val) ...
          - 2 * sigma / y(1) - 4 * mu * y(2) / y(1)) / rho - 1.5 * y(2)^2) / y(1) ...
    ];
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
    [~, y] = ode15s(rp_matlab, t, [R0; 0.0], opts);
end

R     = y(:, 1);          % Bubble radius [m]
R_dot = y(:, 2);          % Wall velocity [m/s]
Mw    = abs(R_dot) / c;   % Bubble-wall Mach number [-]

% 4. Extract Primary Oscillation Cycle & Align Main Collapse to t = 0
win_mask = (t >= 4.0e-6) & (t <= 16.5e-6);
t_w  = t(win_mask);
R_w  = R(win_mask);
Mw_w = Mw(win_mask);

[R_max_m, idx_max] = max(R_w);
coll_search = idx_max : min(length(t_w), idx_max + round(0.4 * length(t_w)));
[R_min_m, rel_min] = min(R_w(coll_search));
idx_min = coll_search(rel_min);

t_shift_us   = t_w(idx_min) * 1e6;
t_aligned_us = t_w * 1e6 - t_shift_us;
R_um         = R_w * 1e6;

reb_search = idx_min : min(length(t_w), idx_min + round(0.3 * length(t_w)));
[R_reb_um, ~] = max(R_um(reb_search));

[~, exp_max_idx] = max(R_exp);

fprintf('\n--- Numerical Results (Primary Cycle) ---\n');
fprintf('Applied Time Shift:      %.4f us\n', t_shift_us);
fprintf('R_max:                   %.2f um (Exp: %.2f um)\n', R_max_m*1e6, R_exp(exp_max_idx));
fprintf('Collapse Duration:       %.2f us (Exp: %.2f us)\n', -t_aligned_us(idx_max), -t_exp(exp_max_idx));
fprintf('First Rebound (R_reb):   %.2f um (Exp: %.2f um)\n', R_reb_um, max(R_exp(t_exp > 0.05)));
fprintf('R_min (report separate): %.4f um\n', R_min_m*1e6);
fprintf('Peak Wall Mach Number:   %.2f\n', max(Mw_w(abs(t_aligned_us) < 1.0)));

% 5. Required Plot 1: Experimental-Numerical Radius-Time Overlay (Task c)
f1 = figure(1);
set(f1, 'Position', [100, 100, 760, 440], ...
        'PaperUnits', 'inches', 'PaperPosition', [0, 0, 7.6, 4.4]);
plot(t_aligned_us, R_um, 'k-', 'LineWidth', 1.6);
hold on;
plot(t_exp, R_exp, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
hold off;
xlim([-8.5, 2.5]);
ylim([0, 21]);
xlabel('Shifted Time, t - t_{collapse} (\mus)');
ylabel('Bubble Radius, R (\mum)');
title('Experimental vs. Incompressible Rayleigh-Plesset Radius-Time Comparison');
legend('Incompressible RP Model', 'Experimental Data (Hoeppe et al., Fig. 4)', 'Location', 'northwest');
grid on;
print(f1, '-dpng', '-r300', 'fig1_radius_overlay.png');

% 6. Plot All 3 Consecutive Forcing Cycles on One Continuous Time Axis (Task b)
f2 = figure(2);
set(f2, 'Position', [150, 150, 820, 440], ...
        'PaperUnits', 'inches', 'PaperPosition', [0, 0, 8.2, 4.4]);

hold on;
colors = {'b', [0.85, 0.33, 0.1], [0.47, 0.67, 0.19]};

for k = 1:3
    c_mask = (t >= (k-1)*T) & (t <= k*T);
    plot(t(c_mask) * 1e6, R(c_mask) * 1e6, 'Color', colors{k}, 'LineWidth', 1.6);
end

% Add vertical dashed dividers at period boundaries t = T and t = 2T
y_top = 60;
line([T*1e6, T*1e6], [0, y_top], 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1.0);
line([2*T*1e6, 2*T*1e6], [0, y_top], 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1.0);
hold off;

xlim([0, 3*T*1e6]);
ylim([0, y_top]);
xlabel('Time, t (\mus)');
ylabel('Bubble Radius, R (\mum)');
title('Bubble Radius Evolution Across Three Consecutive Forcing Cycles');
legend('Forcing Cycle 1 (0 - T)', 'Forcing Cycle 2 (T - 2T)', 'Forcing Cycle 3 (2T - 3T)', 'Location', 'northwest');
grid on;

print(f2, '-dpng', '-r300', 'fig2_continuous_cycles.png');

% 7. Wall Mach Number & Compressibility Breakdown (Task e)
f3 = figure(3);
set(f3, 'Position', [200, 200, 760, 440], ...
        'PaperUnits', 'inches', 'PaperPosition', [0, 0, 7.6, 4.4]);

semilogy(t_aligned_us, Mw_w, 'b-', 'LineWidth', 1.4);
hold on;
line([-8.5, 2.5], [0.3, 0.3], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.2);
hold off;
xlim([-8.5, 2.5]);
ylim([1e-4, 20]);
xlabel('Shifted Time, t - t_{collapse} (\mus)');
ylabel('Wall Mach Number, M_w = |dR/dt|/c');
title('Bubble-Wall Mach Number and Compressibility Breakdown');
legend('Wall Mach Number M_w(t)', 'Compressibility Limit (M_w = 0.3)', 'Location', 'northwest');
grid on;

print(f3, '-dpng', '-r300', 'fig3_wall_mach.png');
fprintf('Done! Saved fig1_radius_overlay.png, fig2_continuous_cycles.png, and fig3_wall_mach.png.\n');
