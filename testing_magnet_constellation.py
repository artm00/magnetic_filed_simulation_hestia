import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from create_magnet_constellation import create_magnet_constellation


def print_constellation(array):
    symbol_map = [' ', '↑', '↓', '←', '→', 'X', 'O']
    rows, cols = array.shape
    for i in range(rows):
        row_str = ''.join(symbol_map[array[i, j]] + ' ' for j in range(cols))
        print(row_str)


def visualize_constellation(array):
    rows, cols = array.shape

    colors = [
        [1.0, 1.0, 1.0],   # 0: empty  - white
        [0.2, 0.8, 0.2],   # 1: up     - green
        [0.8, 0.2, 0.2],   # 2: down   - red
        [0.2, 0.2, 0.8],   # 3: left   - blue
        [0.8, 0.8, 0.2],   # 4: right  - yellow
        [0.5, 0.5, 0.5],   # 5: X      - gray
        [0.9, 0.9, 0.9],   # 6: O      - light gray
    ]
    cmap = mcolors.ListedColormap(colors)
    bounds = np.arange(-0.5, 7.5, 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(facecolor='white')
    im = ax.imshow(array, cmap=cmap, norm=norm, origin='upper')
    ax.set_aspect('equal')

    for i in range(rows + 1):
        ax.axhline(i - 0.5, color='black', linewidth=0.5)
    for j in range(cols + 1):
        ax.axvline(j - 0.5, color='black', linewidth=0.5)

    for i in range(rows):
        for j in range(cols):
            val = array[i, j]
            if val == 1:    # up arrow
                ax.quiver(j, i, 0, -0.3, scale=1, scale_units='xy',
                          angles='xy', color='black', linewidth=2)
            elif val == 2:  # down arrow
                ax.quiver(j, i, 0, 0.3, scale=1, scale_units='xy',
                          angles='xy', color='black', linewidth=2)
            elif val == 3:  # left arrow
                ax.quiver(j, i, -0.3, 0, scale=1, scale_units='xy',
                          angles='xy', color='black', linewidth=2)
            elif val == 4:  # right arrow
                ax.quiver(j, i, 0.3, 0, scale=1, scale_units='xy',
                          angles='xy', color='black', linewidth=2)
            elif val == 5:  # X
                ax.text(j, i, 'X', ha='center', va='center',
                        fontsize=12, fontweight='bold')
            elif val == 6:  # O
                ax.text(j, i, 'O', ha='center', va='center',
                        fontsize=12, fontweight='bold')

    ax.set_title(f'Magnet Constellation Array ({rows} x {cols})')
    ax.set_xlabel('Column')
    ax.set_ylabel('Row')

    cbar = fig.colorbar(im, ax=ax, ticks=range(7))
    cbar.ax.set_yticklabels(['Empty', 'Up', 'Down', 'Left', 'Right', 'X', 'O'])

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    n_z_magnets = 6
    constellation = create_magnet_constellation(n_z_magnets)

    print(f'Array shape: {constellation.shape[0]} x {constellation.shape[1]}\n')

    print('Magnet constellation (numeric):')
    print(constellation)

    print('\nMagnet constellation (symbols):')
    print_constellation(constellation)

    visualize_constellation(constellation)

    np.savetxt('magnet_constellation.txt', constellation, fmt='%d')
    np.save('magnet_constellation.npy', constellation)
