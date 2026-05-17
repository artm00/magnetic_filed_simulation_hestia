# Replication of Xu's Figure 3.5 - Bx Component Only
# First-order harmonic model (Equation 3.7)
#
import numpy as np
import matplotlib.pyplot as plt

from magnetic_nodes_halbach import magnetic_nodes_halbach


# --- System Parameters (from Xu's thesis Section 3.3.3) ---
tau = 50.8e-3              # Pole pitch (m)
tau_m = 25.4e-3            # Pole width (m)
pole_width = 2 * 0.0254    # Pole: 2 inches to metres (overwritten below)
pole_width = tau / 2
hm = 0.5 * 0.0254          # Magnet thickness: 0.5 inches to metres
Br = 1.32                  # Remanence (T), typical NdFeB

# Test line parameters
x_start = 3.5e-3           # Start: 3.5 mm
y_start = 10.5e-3          # Start: 10.5 mm
x_end = 192.5e-3           # End: 192.5 mm
y_end = 37.5e-3            # End: 37.5 mm
z = -25e-3  # Fixed evaluation height

mu0 = 4 * np.pi * 1e-7    # Permeability of free space

# --- Generate Test Line Points ---
num_points = 200
x_line = np.linspace(x_start, x_end, num_points)
y_line = np.linspace(y_start, y_end, num_points)

# --- K' coefficient (Equation 3.7) ---
mt = hm   # Top of magnet array
mb = 0    # Bottom of magnet array at z=0

K_prime = (
    (4 * Br / (np.pi**4 * mu0))
    * (np.exp(-np.sqrt(2) * np.pi * mt / tau) - np.exp(-np.sqrt(2) * np.pi * mb / tau))
    * (np.pi * np.sin(tau_m * np.pi / tau) - np.sqrt(2) * np.cos(tau_m * np.pi / tau) + np.sqrt(2) * tau)
)

lambda_ = np.sqrt(2) * np.pi / tau  # Wave number

# --- Calculate Bx field at each point ---
Bx_harmonic = np.zeros(num_points)
Bx_nodes = np.zeros(num_points)

for i in range(num_points):
    x = x_line[i]
    y = y_line[i]

    decay = np.exp(lambda_ * z)  # z is negative

    Bx_harmonic[i] = -mu0 * K_prime * decay * np.cos(np.pi * x / tau) * np.sin(np.pi * y / tau)
    Bx_nodes[i], _, _ = magnetic_nodes_halbach(x, y, z, tau, tau_m, hm, Br)

# --- Convert x-position to mm ---
x_mm = x_line * 1000

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(x_mm, Bx_harmonic, 'b-', linewidth=2, label='Harmonic (Eq 3.7)')
ax.plot(x_mm, Bx_nodes, 'r--', linewidth=2, label='Magnetic Nodes')
ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
ax.set_xlabel('x (mm)', fontsize=14)
ax.set_ylabel('B_x (T)', fontsize=14)
ax.set_title('X-component of Magnetic Flux Density - First-Order Harmonic', fontsize=14)
ax.grid(True)
ax.set_ylim([-0.175, 0.175])
ax.set_xlim([0, 200])
ax.legend()

plt.tight_layout()
plt.show()

print('\n=== Bx Field Results ===')
print(f'Peak Bx: {np.max(np.abs(Bx_harmonic)):.4f} T')
print(f'Min Bx:  {np.min(Bx_harmonic):.4f} T')
print(f'Max Bx:  {np.max(Bx_harmonic):.4f} T')
