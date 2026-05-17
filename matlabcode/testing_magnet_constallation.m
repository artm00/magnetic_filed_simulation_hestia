% Example script to test the function
clear; clc;

% Create constellation for 6 Z magnets
n_z_magnets = 6;
constellation = create_magnet_constellation(n_z_magnets);

% Display array size
fprintf('Array shape: %d x %d\n\n', size(constellation, 1), size(constellation, 2));

% Display numeric array
disp('Magnet constellation (numeric):');
disp(constellation);

% Display with symbols
fprintf('\nMagnet constellation (symbols):\n');
print_constellation(constellation);

% Visualize the constellation
visualize_constellation(constellation);

% Save to file
writematrix(constellation, 'magnet_constellation.txt');

% You can also save as .mat file
save('magnet_constellation.mat', 'constellation');



function print_constellation(array)
    % PRINT_CONSTELLATION Print the constellation array with symbols for visualization
    
    symbol_map = {' ', '↑', '↓', '←', '→', 'X', 'O'};
    
    [rows, cols] = size(array);
    
    for i = 1:rows
        row_str = '';
        for j = 1:cols
            val = array(i, j) + 1;  % +1 because MATLAB is 1-indexed
            row_str = [row_str, symbol_map{val}, ' '];
        end
        fprintf('%s\n', row_str);
    end
end

function visualize_constellation(array)
    % VISUALIZE_CONSTELLATION Create a visual representation of the magnet array
    
    figure('Color', 'white');
    imagesc(array);
    colormap([1 1 1; 0.2 0.8 0.2; 0.8 0.2 0.2; 0.2 0.2 0.8; 0.8 0.8 0.2; 0.5 0.5 0.5; 0.9 0.9 0.9]);
    axis equal tight;
    
    % Add grid
    hold on;
    [rows, cols] = size(array);
    for i = 0.5:rows+0.5
        plot([0.5, cols+0.5], [i, i], 'k-', 'LineWidth', 0.5);
    end
    for j = 0.5:cols+0.5
        plot([j, j], [0.5, rows+0.5], 'k-', 'LineWidth', 0.5);
    end
    
    % Add arrows and symbols
    for i = 1:rows
        for j = 1:cols
            val = array(i, j);
            if val == 1  % up arrow
                quiver(j, i, 0, -0.3, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 1);
            elseif val == 2  % down arrow
                quiver(j, i, 0, 0.3, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 1);
            elseif val == 3  % left arrow
                quiver(j, i, -0.3, 0, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 1);
            elseif val == 4  % right arrow
                quiver(j, i, 0.3, 0, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 1);
            elseif val == 5  % X
                text(j, i, 'X', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
            elseif val == 6  % O
                text(j, i, 'O', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
            end
        end
    end
    
    title(sprintf('Magnet Constellation Array (%d x %d)', rows, cols));
    xlabel('Column');
    ylabel('Row');
    colorbar('Ticks', [0:6], 'TickLabels', {'Empty', 'Up', 'Down', 'Left', 'Right', 'X', 'O'});
end