function [Bx, By, Bz] = magnetic_nodes_halbach(x, y, z, tau, tau_m, hm, Br)
% MAGNETIC_NODES_HALBACH Calculate magnetic field using magnetic nodes method
%
% [Bx, By, Bz] = magnetic_nodes_halbach(x, y, z, tau, tau_m, hm, Br)
%
% Inputs:
%   x, y, z - Coordinates of point of interest (in meters)
%             Origin at center of bottom surface of Halbach array
%   tau     - Pole pitch (m)
%   tau_m   - Z magnet width (m)
%   hm      - Magnet height (m)
%   Br      - Remanence (T)
%
% Outputs:
%   Bx, By, Bz - Magnetic flux density components (in Tesla)
%
% Based on Xu's thesis Equations 3.3-3.5 and magnetic nodes method
% Models a 2D Halbach array using create_magnet_constellation pattern
    
    % For a 6x6 Z-magnet array, we need n_z_magnets = 3
    % This creates a 2*3+1 = 7x7 constellation (with some empty positions)
    % But Xu has 6x6 actual magnets, so we'll use n_z_magnets = 3
    n_z_magnets = 3;
    
    %% Generate Halbach Array Layout using constellation function
    constellation = create_magnet_constellation(n_z_magnets);
    array_size = size(constellation, 1);  % Should be 2*n_z_magnets + 1
    
    % Initialize field components
    Bx_total = 0;
    By_total = 0;
    Bz_total = 0;
    
    % Calculate magnet spacing
    % For a Halbach array with pole pitch tau, magnets are spaced at tau/2
    magnet_spacing = tau / 2;
    
    % Loop through all positions in constellation
    for i = 1:array_size
        for j = 1:array_size
            % Get magnet type at this position
            magnet_type = constellation(i, j);
            
            % Skip if no magnet (0 = empty)
            if magnet_type == 0
                continue;
            end
            
            % Calculate center position of this magnet
            % Array is centered at origin
            magnet_center_x = (i - (array_size+1)/2) * magnet_spacing;
            magnet_center_y = (j - (array_size+1)/2) * magnet_spacing;
            magnet_center_z = hm / 2;  % Center height of magnet
            
            % Map constellation codes to magnetization directions
            mag_dir = map_constellation_to_direction(magnet_type);
            
            % Determine magnet dimensions based on type
            % Z-magnets (up/down) use tau_m width
            % X/Y magnets use smaller dimensions
            if magnet_type == 1 || magnet_type == 2  % Up or Down (Z-magnets)
                wx = tau_m;
                wy = tau_m;
            elseif magnet_type == 3 || magnet_type == 4  % Left or Right (X-magnets)
                wx = (tau - tau_m) / 2;  % Transitional magnet width
                wy = tau_m;
            elseif magnet_type == 5 || magnet_type == 6  % X or O (diagonal magnets)
                wx = (tau - tau_m) / 2;
                wy = (tau - tau_m) / 2;
            else
                continue;  % Skip unknown types
            end
            
            % Get the 8 corner nodes and their contributions
            [Bx_mag, By_mag, Bz_mag] = calculate_magnet_field(x, y, z, ...
                magnet_center_x, magnet_center_y, magnet_center_z, ...
                wx, wy, hm, Br, mag_dir);
            
            % Accumulate contributions
            Bx_total = Bx_total + Bx_mag;
            By_total = By_total + By_mag;
            Bz_total = Bz_total + Bz_mag;
        end
    end
    
    % Return total field
    Bx = Bx_total;
    By = By_total;
    Bz = Bz_total;
end

%% Helper function: Map constellation code to magnetization direction
function mag_dir = map_constellation_to_direction(code)
% Map constellation codes to magnetization directions
% 1 = up (↑)    → Z+
% 2 = down (↓)  → Z-
% 3 = left (←)  → X-
% 4 = right (→) → X+
% 5 = X (crossed) → diagonal NE-SW
% 6 = O (circle)  → diagonal NW-SE

    switch code
        case 1
            mag_dir = 'Z+';  % Up
        case 2
            mag_dir = 'Z-';  % Down
        case 3
            mag_dir = 'X-';  % Left
        case 4
            mag_dir = 'X+';  % Right
        case 5
            mag_dir = 'XY+';  % X position - diagonal (+X+Y)
        case 6
            mag_dir = 'XY-';  % O position - diagonal (+X-Y) or (-X+Y)
        otherwise
            mag_dir = 'Z+';  % Default
    end
end

%% Helper function: Calculate field from a single magnet using 8 nodes
function [Bx, By, Bz] = calculate_magnet_field(xm, ym, zm, ...
    cx, cy, cz, wx, wy, hz, Br, mag_dir)
% Calculate magnetic field from a single cuboidal magnet
% (xm, ym, zm) - measurement point
% (cx, cy, cz) - magnet center
% (wx, wy, hz) - magnet width in x, y, and height in z
% Br - remanence
% mag_dir - magnetization direction: 'X+', 'X-', 'Y+', 'Y-', 'Z+', 'Z-', 'XY+', 'XY-'

    % Define the 8 corner positions relative to magnet center
    corners = [
        cx - wx/2, cy - wy/2, cz - hz/2;  % 1: bottom-front-left
        cx + wx/2, cy - wy/2, cz - hz/2;  % 2: bottom-front-right
        cx - wx/2, cy + wy/2, cz - hz/2;  % 3: bottom-back-left
        cx + wx/2, cy + wy/2, cz - hz/2;  % 4: bottom-back-right
        cx - wx/2, cy - wy/2, cz + hz/2;  % 5: top-front-left
        cx + wx/2, cy - wy/2, cz + hz/2;  % 6: top-front-right
        cx - wx/2, cy + wy/2, cz + hz/2;  % 7: top-back-left
        cx + wx/2, cy + wy/2, cz + hz/2;  % 8: top-back-right
    ];
    
    % Determine node charges based on magnetization direction
    epsilon = zeros(8, 1);
    
    switch mag_dir
        case 'X+'  % Magnetized in +X direction
            epsilon([1,3,5,7]) = -Br;  % Left face (x = cx - wx/2)
            epsilon([2,4,6,8]) = +Br;  % Right face (x = cx + wx/2)
        case 'X-'  % Magnetized in -X direction
            epsilon([1,3,5,7]) = +Br;
            epsilon([2,4,6,8]) = -Br;
        case 'Y+'  % Magnetized in +Y direction
            epsilon([1,2,5,6]) = -Br;  % Front face (y = cy - wy/2)
            epsilon([3,4,7,8]) = +Br;  % Back face (y = cy + wy/2)
        case 'Y-'  % Magnetized in -Y direction
            epsilon([1,2,5,6]) = +Br;
            epsilon([3,4,7,8]) = -Br;
        case 'Z+'  % Magnetized in +Z direction
            epsilon([1,2,3,4]) = -Br;  % Bottom face (z = cz - hz/2)
            epsilon([5,6,7,8]) = +Br;  % Top face (z = cz + hz/2)
        case 'Z-'  % Magnetized in -Z direction
            epsilon([1,2,3,4]) = +Br;
            epsilon([5,6,7,8]) = -Br;
        case 'XY+'  % Diagonal +X+Y (northeast)
            % Approximate with combined charges
            epsilon([1,3,5,7]) = -Br/sqrt(2);  % -X face
            epsilon([2,4,6,8]) = +Br/sqrt(2);  % +X face
            epsilon([1,2,5,6]) = epsilon([1,2,5,6]) - Br/sqrt(2);  % -Y face
            epsilon([3,4,7,8]) = epsilon([3,4,7,8]) + Br/sqrt(2);  % +Y face
        case 'XY-'  % Diagonal +X-Y (southeast) or -X+Y (northwest)
            epsilon([1,3,5,7]) = -Br/sqrt(2);  % -X face
            epsilon([2,4,6,8]) = +Br/sqrt(2);  % +X face
            epsilon([1,2,5,6]) = epsilon([1,2,5,6]) + Br/sqrt(2);  % +Y face
            epsilon([3,4,7,8]) = epsilon([3,4,7,8]) - Br/sqrt(2);  % -Y face
    end
    
    % Calculate field contribution from each node (Equations 3.3-3.5)
    Bx = 0;
    By = 0;
    Bz = 0;
    
    for k = 1:8
        xk = corners(k, 1);
        yk = corners(k, 2);
        zk = corners(k, 3);
        
        % Relative distance
        xr = xm - xk;
        yr = ym - yk;
        zr = zm - zk;
        
        % Distance
        r = sqrt(xr^2 + yr^2 + zr^2);
        
        % Avoid singularity at r = 0
        if r < 1e-10
            continue;
        end
        
        % Magnetic node charge
        eps_k = epsilon(k);
        
        % Calculate field components (Equations 3.3, 3.4, 3.5)
        Bxk = (eps_k / (4*pi)) * log(-yr + r);
        Byk = (eps_k / (4*pi)) * log(-xr + r);
        Bzk = (eps_k / (4*pi)) * atan2(xr*yr, zr*r);
        
        % Accumulate
        Bx = Bx + Bxk;
        By = By + Byk;
        Bz = Bz + Bzk;
    end
end