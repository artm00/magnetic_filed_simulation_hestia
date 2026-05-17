function array = create_magnet_constellation(n_z_magnets)
    % CREATE_MAGNET_CONSTELLATION Create a magnet constellation array for a planar motor stator
    %
    % Parameters:
    %   n_z_magnets - Number of magnets in the Z direction (must be even)
    %
    % Returns:
    %   array - Matrix representing magnet orientations where:
    %           1 = up (↑)
    %           2 = down (↓)
    %           3 = left (←)
    %           4 = right (→)
    %           5 = X position (no magnet or crossed)
    %           6 = O position (circular/null)
    %           0 = empty/no magnet
    
    if nargin < 1
        n_z_magnets = 6;  % Default value
    end
    
    % Array size is 2*n_z_magnets + 1
    array_size = 2 * n_z_magnets + 1;
    array = zeros(array_size, array_size);
    
    % Define the repeating 2x2 unit cell pattern
    % This pattern repeats across the array
    
    for i = 1:array_size
        for j = 1:array_size
            % Determine position in 4x4 repeating pattern
            row_mod = mod(i-1, 4);
            col_mod = mod(j-1, 4);
            
            % Pattern assignment based on position
            if row_mod == 0
                if col_mod == 0
                    array(i, j) = 0;  % corner/edge
                elseif col_mod == 1
                    array(i, j) = 1;  % up
                elseif col_mod == 2
                    array(i, j) = 0;  % edge
                elseif col_mod == 3
                    array(i, j) = 2;  % down
                end
            elseif row_mod == 1
                if col_mod == 0
                    array(i, j) = 3;  % left
                elseif col_mod == 1
                    array(i, j) = 5;  % X
                elseif col_mod == 2
                    array(i, j) = 4;  % right
                elseif col_mod == 3
                    array(i, j) = 6;  % O
                end
            elseif row_mod == 2
                if col_mod == 0
                    array(i, j) = 0;  % edge
                elseif col_mod == 1
                    array(i, j) = 2;  % down
                elseif col_mod == 2
                    array(i, j) = 0;  % edge
                elseif col_mod == 3
                    array(i, j) = 1;  % up
                end
            elseif row_mod == 3
                if col_mod == 0
                    array(i, j) = 4;  % right
                elseif col_mod == 1
                    array(i, j) = 6;  % O
                elseif col_mod == 2
                    array(i, j) = 3;  % left
                elseif col_mod == 3
                    array(i, j) = 5;  % X
                end
            end
        end
    end
end

