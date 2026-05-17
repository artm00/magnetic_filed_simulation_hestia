import numpy as np

def compute_b_field(coords, node_positions, epsilon):
    """
    Calculates the 3D magnetic B-field vector at a specific evaluation point.
    Result of the function is magnetic field vector at the point coords in the Magnet coordinate system.
    Parameters:
    -----------
    coords : array-like of shape (3,)
        The target evaluation coordinates [x, y, z] in meters.
    node_positions : array-like of shape (N, 4)
        An array of N source points where each row is [x_k, y_k, z_k, k].
    epsilon : float
        
    Returns:
    --------
    B_vector : ndarray of shape (3,)
        The resulting magnetic field vector [Bx, By, Bz] at the target point.
    """
    k = node_positions[:, 3]
    R = coords - node_positions[:, :3]  # 8 displacement vectors (x, y, z)

    Bx = k * epsilon / (4 * np.pi) * np.log(-R[:,1] + np.sqrt(R[:,0]**2 + R[:,1]**2 + R[:,2]**2))
    By = k * epsilon / (4 * np.pi) * np.log(-R[:,0] + np.sqrt(R[:,0]**2 + R[:,1]**2 + R[:,2]**2))
    Bz = k * epsilon / (4 * np.pi) * np.arctan((R[:,0]*R[:,1]) / (R[:,2]*np.sqrt(R[:,0]**2 + R[:,1]**2 + R[:,2]**2))) ### !!!! maybe z missing, check in the literature again

    B_contributions = np.array([Bx, By, Bz])
    B_vector = np.sum(B_contributions, axis=1)

    return B_vector