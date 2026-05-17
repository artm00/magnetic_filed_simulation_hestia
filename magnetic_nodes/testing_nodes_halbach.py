import numpy as np
from b_field_in_magnet_coords import compute_b_field

# --- System Parameters ---
tau = 50.8e-3
tau_m = 25.4e-3
hm = 0.5 * 0.0254
Br = 1.32

length = tau_m
width = tau_m
height = hm
sign = -1

# 8 corner nodes
pt_0 = np.array([0,      0,     0,      -sign])
pt_1 = np.array([length, 0,     0,      sign])
pt_2 = np.array([0,      width, 0,      sign])
pt_3 = np.array([length, width, 0,      -sign])
pt_4 = np.array([0,      0,     height, sign])
pt_5 = np.array([length, 0,     height, -sign])
pt_6 = np.array([0,      width, height, -sign])
pt_7 = np.array([length, width, height, sign])

node_positions = np.array([pt_0, pt_1, pt_2, pt_3, pt_4, pt_5, pt_6, pt_7])

# Full 4×4 transformation matrix
T = np.array([
    [1,  0,  0, 0],
    [0,  0,  1, 0],
    [0, -1,  0, length],
    [0,  0,  0, 1]
])

# T = np.array([
#     [1,  0,  0, 0],
#     [0,  1,  0, 0],
#     [0, 0,  1, 0],
#     [0,  0,  0, 1]
# ])

# Extract rotation matrix (top-left 3×3) and its inverse
R = T[:3, :3]  # rotation part only
R_inv = np.linalg.inv(R)  # inverse rotation

# Extract translation (for transforming positions)
translation = T[:3, 3]

# ------------------------- Computing B field -------------------------
nx, ny, nz = 10, 10, 7
x_coords = np.linspace(-2 * length, 3 * length, nx)
y_coords = np.linspace(-2 * width,  3 * width,  ny)
z_coords = np.linspace(-height,     2 * height, nz)

XXX, YYY, ZZZ = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')

B_volume = np.zeros((nx, ny, nz, 3))
for i, x in enumerate(x_coords):
    for j, y in enumerate(y_coords):
        for k, z in enumerate(z_coords):
            # Position in world coordinates
            coord_world = np.array([x, y, z])
            
            # Transform position to magnet coordinates
            # coord_mag = R_inv @ (coord_world - translation)
            coord_mag = R_inv @ coord_world - R_inv @ translation
            
            # Compute B field in magnet coordinate system
            B_mag = compute_b_field(coords=coord_mag, 
                                   node_positions=node_positions, 
                                   epsilon=Br)
            B_mag = np.array(B_mag)
            
            # Transform B field vector back to world coordinates
            # Only apply rotation, NOT translation
            B_global = R @ B_mag
            
            B_volume[i, j, k] = B_global

# ------------------------- PLOTTING -------------------------
import matplotlib.pyplot as plt

U = B_volume[:, :, :, 0]
V = B_volume[:, :, :, 1]
W = B_volume[:, :, :, 2]

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection='3d')

ax.quiver(XXX * 1000, YYY * 1000, ZZZ * 1000,
          U, V, W,
          length=4, normalize=True, linewidth=0.6, alpha=0.7)

# Transform magnet corners to world coordinates for plotting
corners_mag = np.array([[0,0,0],[length,0,0],[0,width,0],[length,width,0],
                        [0,0,height],[length,0,height],[0,width,height],[length,width,height]])
corners_world = (R @ corners_mag.T).T + translation
corners_world *= 1000  # to mm

edges = [(0,1),(0,2),(1,3),(2,3),(4,5),(4,6),(5,7),(6,7),(0,4),(1,5),(2,6),(3,7)]
for a, b in edges:
    ax.plot(*zip(corners_world[a], corners_world[b]), color='red', linewidth=1.5)

ax.set_xlabel('x (mm)')
ax.set_ylabel('y (mm)')
ax.set_zlabel('z (mm)')
ax.set_title('B field in 3D (Bx, By, Bz arrows)')
plt.tight_layout()
plt.show()