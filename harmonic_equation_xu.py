# Replication of Xu's Figure 3.5 - Harmonic Model (Equation 3.7)
# Comparison of first-order harmonic magnetic field prediction
#
# NOTE: This script is an incomplete translation of harmonic_equasion_Xu.m.
# In the original, variables `z` and `y` (used in the field calculation loop)
# were undefined because y_line and z_line were commented out.
# See chatty_version.py for the more complete version.

import numpy as np
import matplotlib.pyplot as plt

# --- System Parameters (from Xu's thesis Section 3.3.3) ---
tau = 50.8e-3        # Pole pitch (m)
tau_m = 25.4e-3      # Pole width (m)
hm = 6.35e-3         # Magnet height (m)
Br = 1.32            # Remanence (T), typical NdFeB

# Test line parameters
x_start = 0
x_end = 192.5e-3     # End: 192.5 mm

mu0 = 4 * np.pi * 1e-7  # Permeability of free space

# --- Generate Test Line Points ---
num_points = 200
x_line = np.linspace(x_start, x_end, num_points)

# --- K' coefficient (Equation 3.7) ---
mt = hm   # Top of magnet array
mb = 0    # Bottom of magnet array at z=0

K_prime = (
    (4 * Br / (np.pi**4 * mu0))
    * (np.exp(-np.sqrt(2) * np.pi * mt / tau) - np.exp(-np.sqrt(2) * np.pi * mb / tau))
    * (np.pi * np.sin(tau_m * np.pi / tau) - np.sqrt(2) * np.cos(tau_m * np.pi / tau) + np.sqrt(2) * tau)
)

# --- Field calculation ---
Bx_harmonic = np.zeros(num_points)

# NOTE: `z` and `y` are not defined in the original script (commented out).
# Assign placeholder values to allow the loop to run.
z = -25e-3  # placeholder: height below magnet
y = 0.0     # placeholder: y-coordinate

for i in range(num_points):
    x = x_line[i]

    decay = np.exp(np.sqrt(2) * np.pi / tau * z)  # z is negative

    Bx_harmonic[i] = -mu0 * K_prime * decay * np.cos(np.pi * x / tau) * np.sin(np.pi * y / tau)

# --- Plot ---
fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(x_line * 1000, Bx_harmonic, 'b-', linewidth=2, label='First-order Harmonic (Eq 3.7)')
ax.set_xlabel('Distance along line (mm)', fontsize=12)
ax.set_ylabel('B_x (T)', fontsize=12)
ax.set_title('X-component of Magnetic Flux Density', fontsize=13)
ax.grid(True)
ax.set_ylim([-0.175, 0.175])
ax.legend(loc='best')

plt.tight_layout()
plt.show()
