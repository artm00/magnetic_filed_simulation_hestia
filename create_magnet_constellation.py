import numpy as np


def create_magnet_constellation(n_z_magnets=6):
    """
    Create a magnet constellation array for a planar motor stator.

    Parameters:
        n_z_magnets: Number of magnets in the Z direction (must be even)

    Returns:
        array: 2D numpy array with magnet orientations:
               0=empty, 1=up(↑), 2=down(↓), 3=left(←), 4=right(→), 5=X, 6=O
    """
    array_size = 2 * n_z_magnets + 1
    array = np.zeros((array_size, array_size), dtype=int)

    pattern = {
        (0, 0): 0, (0, 1): 1, (0, 2): 0, (0, 3): 2,
        (1, 0): 3, (1, 1): 5, (1, 2): 4, (1, 3): 6,
        (2, 0): 0, (2, 1): 2, (2, 2): 0, (2, 3): 1,
        (3, 0): 4, (3, 1): 6, (3, 2): 3, (3, 3): 5,
    }

    for i in range(array_size):
        for j in range(array_size):
            array[i, j] = pattern[(i % 4, j % 4)]

    return array
