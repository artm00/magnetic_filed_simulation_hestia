import numpy as np


def compute_b_harmonic_first_order(points, tau, tau_m, hm, Br):
    """
    Calculate B field using first-order harmonic model (Xu Equation 3.7).
    
    Parameters:
    -----------
    points : ndarray (N, 3)
        Sample points in world coordinates [x, y, z]
    tau : float
        Pole pitch (m)
    tau_m : float
        Pole width (m)
    hm : float
        Magnet thickness (m)
    Br : float
        Remanence (T)
    
    Returns:
    --------
    B_harmonic : ndarray (N, 3)
        Magnetic field [Bx, By, Bz] at each point
    """
    mu0 = 4 * np.pi * 1e-7  # Permeability of free space
    
    # Calculate K' coefficient (Equation 3.7)
    mt = hm  # Top of magnet array
    mb = 0   # Bottom of magnet array at z=0
    
    K_prime = (4*Br / (np.pi**4 * mu0)) * \
              (np.exp(-np.sqrt(2)*np.pi*mt/tau) - np.exp(-np.sqrt(2)*np.pi*mb/tau)) * \
              (np.pi * np.sin(tau_m*np.pi/tau) - np.sqrt(2)*np.cos(tau_m*np.pi/tau) + np.sqrt(2)*tau)
    
    # Wave number
    lam = np.sqrt(2) * np.pi / tau
    
    # Initialize output
    B_harmonic = np.zeros((len(points), 3))
    
    for i, pt in enumerate(points):
        x, y, z = pt
        
        # Exponential decay with height (z is typically negative)
        decay = np.exp(lam * z)
        
        # Equation 3.7 components
        Bx = -mu0 * K_prime * decay * np.cos(np.pi*x/tau) * np.sin(np.pi*y/tau)
        By = -mu0 * K_prime * decay * np.sin(np.pi*x/tau) * np.cos(np.pi*y/tau)
        Bz = -mu0 * K_prime * decay * np.sqrt(2) * np.sin(np.pi*x/tau) * np.sin(np.pi*y/tau)
        
        B_harmonic[i] = [Bx, By, Bz]
    
    return B_harmonic