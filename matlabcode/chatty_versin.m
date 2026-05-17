%% Replication of Xu's Figure 3.5 - Bx Component Only
% First-order harmonic model (Equation 3.7)
clear; close all; clc;

%% System Parameters (from Xu's thesis Section 3.3.3)
% 6x6 Halbach array specifications

tau = 50.8e-3;           % Pole pitch: 
tau_m = 25.4e-3;      % Pole:
pole_width = 2 * 0.0254;      % Pole: 2 inches to meters  
pole_width = tau / 2;
hm = 0.5 * 0.0254;            % Magnet thickness: 0.5 inches to meters
Br = 1.32;                    % Remanence: 1.32 T (typical NdFeB)

% Test line parameters
x_start = 3.5e-3;             % Start: 3.5 mm
y_start = 10.5e-3;            % Start: 10.5 mm  
x_end = 192.5e-3;             % End: 192.5 mm
y_end = 37.5e-3;              % End: 37.5 mm
z_height = -25e-3;            % Height: -25 mm below magnet

% Physical constants
mu0 = 4*pi*1e-7;              % Permeability of free space

%% Generate Test Line Points
num_points = 200;
x_line = linspace(x_start, x_end, num_points);
y_line = linspace(y_start, y_end, num_points);
z_line = z_height * ones(1, num_points);

%% Calculate K' coefficient (Equation 3.7)
mt = hm;             % Top of magnet array
mb = 0;              % Bottom of magnet array at z=0

K_prime = (4*Br / (pi^4 * mu0)) * ...
          (exp(-sqrt(2)*pi*mt/tau) - exp(-sqrt(2)*pi*mb/tau)) * ...
          (pi * sin(tau_m*pi/tau) - sqrt(2)*cos(tau_m*pi/tau) + sqrt(2)*tau);

% Wave number
lambda = sqrt(2) * pi / tau;

%% Calculate Bx field at each point
Bx_harmonic = zeros(1, num_points);
Bx_nodes = zeros(1, num_points);


z = -5e-3

for i = 1:num_points
    x = x_line(i);
    y = y_line(i);
    % z = z_line(i);
    
    
    % Exponential decay with height
    decay = exp(lambda * z);  % z is negative
    
    % Bx component from Equation 3.7
    Bx_harmonic(i) = -mu0 * K_prime * decay * cos(pi*x/tau) * sin(pi*y/tau);
    [Bx_nodes(i), ~, ~] = magnetic_nodes_halbach(x, y, z, tau, tau_m, hm, Br);
end

%% Convert x-position to mm for plotting
x_mm = x_line * 1000;  % Convert to mm

%% Plot Results
figure('Position', [100, 100, 800, 500]);
hold on;
plot(x_mm, Bx_harmonic, 'b-', 'LineWidth', 2, 'DisplayName', 'Harmonic (Eq 3.7)');
plot(x_mm, Bx_nodes, 'r--', 'LineWidth', 2, 'DisplayName', 'Magnetic Nodes');
xlabel('x (mm)', 'FontSize', 14);
ylabel('B_x (T)', 'FontSize', 14);
title('X-component of Magnetic Flux Density - First-Order Harmonic', 'FontSize', 14);
grid on;
ylim([-0.175, 0.175]);
xlim([0, 200]);

% Add reference lines
hold on;
plot([0 200], [0 0], 'k--', 'LineWidth', 0.5);

% Print statistics
fprintf('\n=== Bx Field Results ===\n');
fprintf('Peak Bx: %.4f T\n', max(abs(Bx_harmonic)));
fprintf('Min Bx:  %.4f T\n', min(Bx_harmonic));
fprintf('Max Bx:  %.4f T\n', max(Bx_harmonic));