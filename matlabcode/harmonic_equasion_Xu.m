%% Replication of Xu's Figure 3.5 - Harmonic Model (Equation 3.7)
% Comparison of first-order harmonic magnetic field prediction
clear; close all; clc;

%% System Parameters (from Xu's thesis Section 3.3.3)
% 6x6 Halbach array specifications
tau = 50.8e-3;           % Pole pitch: 
tau_m = 25.4e-3;      % Pole:  
hm = 6.35e-3;            % Magnet height 

Br = 1.32;                    % Remanence: 1.32 T (typical NdFeB)

% Test line parameters
x_start = 0;             % Start: 3.5 mm
% y_start = 10.5e-3;            % Start: 10.5 mm  
x_end = 192.5e-3;             % End: 192.5 mm
% y_end = 37.5e-3;              % End: 37.5 mm
% z_height = -25e-3;            % Height: -25 mm below magnet

% Physical constants
mu0 = 4*pi*1e-7;              % Permeability of free space

%% Generate Test Line Points
num_points = 200;
x_line = linspace(x_start, x_end, num_points);
% y_line = linspace(y_start, y_end, num_points);
% z_line = z_height * ones(1, num_points);

%% Calculate Magnetic Field using Harmonic Model (Equation 3.7)
% This is the FIRST-ORDER harmonic approximation

% Calculate K' coefficient (Equation 3.7)
mt = hm;             % Top of magnet array
mb = 0;              % Bottom of magnet array at z=0

% K' calculation - note the complex expression from Eq 3.7
K_prime = (4*Br / (pi^4 * mu0)) * ...
          (exp(-sqrt(2)*pi*mt/tau) - exp(-sqrt(2)*pi*mb/tau)) * ...
          (pi * sin(tau_m*pi/tau) - sqrt(2)*cos(tau_m*pi/tau) + sqrt(2)*tau);

% Initialize field arrays
Bx_harmonic = zeros(1, num_points);
% By_harmonic = zeros(1, num_points);
% Bz_harmonic = zeros(1, num_points);

% Calculate field at each point along the line
for i = 1:num_points
    x = x_line(i);
    % y = y_line(i);
    % z = z_line(i);
    
    % Exponential decay with height
    decay = exp(sqrt(2) * pi / tau * z);  % z is negative
    
    % Harmonic model (Equation 3.7)
    Bx_harmonic(i) = -mu0 * K_prime * decay * cos(pi*x/tau) * sin(pi*y/tau);
    % By_harmonic(i) = -mu0 * K_prime * decay * sin(pi*x/tau) * cos(pi*y/tau);
    % Bz_harmonic(i) = -mu0 * K_prime * decay * sqrt(2) * sin(pi*x/tau) * cos(pi*y/tau);
end









%% Convert distance along line to mm for plotting
distance_along_line = sqrt((x_line - x_start).^2 + (y_line - y_start).^2) * 1000;  % in mm

%% Plot Results - Replicating Figure 3.5 style
figure('Position', [100, 100, 1200, 400]);

% Plot Bx component
subplot(1,2,1);
plot(distance_along_line, Bx_harmonic, 'b-', 'LineWidth', 2);
hold on;
xlabel('Distance along line (mm)', 'FontSize', 12);
ylabel('B_x (T)', 'FontSize', 12);
title('X-component of Magnetic Flux Density', 'FontSize', 13);
grid on;
ylim([-0.175, 0.175]);
legend('First-order Harmonic (Eq 3.7)', 'Location', 'best');

% % Plot Bz component  
% subplot(1,2,2);
% plot(distance_along_line, Bz_harmonic, 'r-', 'LineWidth', 2);
% hold on;
% xlabel('Distance along line (mm)', 'FontSize', 12);
% ylabel('B_z (T)', 'FontSize', 12);
% title('Z-component of Magnetic Flux Density', 'FontSize', 13);
% grid on;
% ylim([-0.8, 0.8]);
% legend('First-order Harmonic (Eq 3.7)', 'Location', 'best');
% 
% sgtitle('Replication of Xu Figure 3.5 - First-Order Harmonic Model', 'FontSize', 14, 'FontWeight', 'bold');

% %% Print some statistics
% fprintf('\n=== Harmonic Model Results ===\n');
% fprintf('Peak Bx: %.4f T\n', max(abs(Bx_harmonic)));
% fprintf('Peak Bz: %.4f T\n', max(abs(Bz_harmonic)));
% fprintf('Central region (0-120mm) captured: %s\n', 'Yes');
% fprintf('\nNote: Harmonic model is only valid in central region.\n');
% fprintf('Edge effects (near 0mm and >150mm) will show deviation from reality.\n');
% 
% %% Optional: Plot field magnitude
% figure('Position', [100, 600, 800, 400]);
% B_magnitude = sqrt(Bx_harmonic.^2 + By_harmonic.^2 + Bz_harmonic.^2);
% plot(distance_along_line, B_magnitude, 'k-', 'LineWidth', 2);
% xlabel('Distance along line (mm)', 'FontSize', 12);
% ylabel('|B| (T)', 'FontSize', 12);
% title('Magnetic Field Magnitude along Test Line', 'FontSize', 13);
% grid on;
% legend('First-order Harmonic');
% 
% %% Display key parameters used
% fprintf('\n=== Parameters Used ===\n');
% fprintf('Pole pitch (tau): %.1f mm\n', tau*1000);
% fprintf('Magnet thickness (hm): %.1f mm\n', hm*1000);
% fprintf('Remanence (Br): %.2f T\n', Br);
% fprintf('Test height: %.1f mm below magnet\n', abs(z_height)*1000);
% fprintf('K_prime coefficient: %.6e\n', K_prime);
% %fprintf('Lambda (decay constant): %.2f m^-1\n', lambda);