import numpy as np
from b_field_in_magnet_coords import compute_b_field
import sys
from pathlib import Path

# Get the absolute path of the directory one level up
parent_dir = str(Path(__file__).resolve().parent.parent)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)
from create_magnet_constellation import create_magnet_constellation






# --- System Parameters (from Xu's thesis Section 3.3.3) ---
tau = 50.8e-3              # Pole pitch (m) = length z magnet plus width x/y magnet
tau_m = 25.4e-3            # Pole width (m) = length/width of z manget
hm = 0.5 * 0.0254          # Magnet thickness: 0.5 inches to metres
Br = 1.32                  # Remanence (T), typical NdFeB

NR_Z_MAGNETS = 1

length = tau_m
width = tau_m
height = hm

# 8 corner nodes of the cuboid; pt_0 is at the origin
pt_0 = np.array([0,      0,     0      , -1])
pt_1 = np.array([length, 0,     0      , 1])
pt_2 = np.array([0,      width, 0      , 1])
pt_3 = np.array([length, width, 0      , -1])
pt_4 = np.array([0,      0,     height , 1])
pt_5 = np.array([length, 0,     height , -1])
pt_6 = np.array([0,      width, height , -1])
pt_7 = np.array([length, width, height , 1])

node_positions = np.array([pt_0, pt_1, pt_2, pt_3, pt_4, pt_5, pt_6, pt_7])

# -------------------------  -------------------------
# Create Trafo matrix for the magnet coordinate sytstem  
# T = np.array([
#     [1,  0,  0, 0],
#     [0,  0,  1, 0],
#     [0, -1,  0, length],
#     [0,  0,  0, 1]
# ])

Trafo_dict = {
    "R_Mx": np.array([
    [0.0,  0.0,  1.0, 0.0],
    [1.0,  0.0,  0.0, 0.0],
    [0.0, 1.0,  0.0, 0.0],
    [0.0,  0.0,  0.0, 1.0]]),

    "R_My": np.array([
    [0.0,  1.0,  0.0, 0.0],
    [0.0,  0.0,  1.0, 0.0],
    [1.0, 0.0,  0.0, 0.0],
    [0.0,  0.0,  0.0, 1.0]]),

    "R_Mz": np.array([
    [1.0,  0.0,  0.0, 0.0],
    [0.0,  1.0,  0.0, 0.0],
    [0.0, 0.0,  1.0, 0.0],
    [0.0,  0.0,  0.0, 1.0]])
}

# ------------------------- Computing B field -------------------------

# constellation first — its size drives the grid extents
constallation = create_magnet_constellation(NR_Z_MAGNETS)
print(constallation)
rows, cols = constallation.shape

# 3D grid covering the full constellation footprint + one-cell buffer on each side
nx, ny, nz = 10, 10, 7
x_coords = np.linspace(-width,        cols * width  + width,  nx)
y_coords = np.linspace(-length,       rows * length + length, ny)
z_coords = np.linspace(-height,       2 * height,             nz)

XXX, YYY, ZZZ = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
B_volume = np.zeros((nx, ny, nz, 3))
for row in range(rows):
    for col in range(cols):
        print(constallation[row, col])
        # selecting the magnetisation direction for the current cell
        sign = 1
        if constallation[row, col] == 0:
            continue  # skip empty cells
        if constallation[row, col] == 1 or constallation[row, col] == 2:
            if constallation[row, col] == 1:
                sign = -1
            T = Trafo_dict["R_My"].copy()
        elif constallation[row, col] == 3 or constallation[row, col] == 4:
            if constallation[row, col] == 3:
                sign = -1
            T = Trafo_dict["R_Mx"].copy()
        elif constallation[row, col] == 5 or constallation[row, col] == 6:
            if constallation[row, col] == 5:
                sign = -1
            T = Trafo_dict["R_Mz"].copy()

        # Creating Translation for the current magnet
        translation = np.array([width*col, length*row, 0, 1])
        T[:, 3] = translation
        T_inv = np.linalg.inv(T)
        Rotation = T[:3, :3]

        for i, x in enumerate(x_coords):
            for j, y in enumerate(y_coords):
                for k, z in enumerate(z_coords):
                    coord = np.array([x, y, z])
                    coord = np.append(coord, 1)
                    # compute coord in the magnet coordinate system
                    
                    coord_mag = T_inv @ coord
                    coord_mag = coord_mag[:3]  # drop homogeneous coordinate

                    # compute B field in the magenet coordinate system
                    B_mag = compute_b_field(coords=coord_mag, node_positions=node_positions, epsilon=sign * Br)
                    B_mag = np.array(B_mag)

                    # transform B field back to global coordinate system
                    B_global = Rotation @ B_mag
                    

                    # if constallation[row, col] == 4:
                    #     B_volume[i, j, k] += B_global
                    B_volume[i, j, k] += B_global
                        

# ------------------------- PLOTTING -------------------------
import matplotlib.pyplot as plt

U = B_volume[:, :, :, 0]
V = B_volume[:, :, :, 1]
W = B_volume[:, :, :, 2]

magnitude = np.sqrt(U**2 + V**2 + W**2).ravel()
norm = plt.Normalize(vmin=magnitude.min(), vmax=magnitude.max())
cmap = plt.cm.plasma
colors = cmap(norm(magnitude))          # (n_arrows, 4) RGBA per arrow

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection='3d')

ax.quiver(XXX * 1000, YYY * 1000, ZZZ * 1000,
          U, V, W,
          length=4, normalize=True, linewidth=0.6, alpha=0.8,
          color=colors)

# colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
fig.colorbar(sm, ax=ax, label='|B| (T)', shrink=0.55, pad=0.1)

# draw each magnet — same loop / T construction as the B field computation
local_corners_hom = np.array([
    [0,      0,     0,      1],
    [length, 0,     0,      1],
    [0,      width, 0,      1],
    [length, width, 0,      1],
    [0,      0,     height, 1],
    [length, 0,     height, 1],
    [0,      width, height, 1],
    [length, width, height, 1],
])
edges = [(0,1),(0,2),(1,3),(2,3),(4,5),(4,6),(5,7),(6,7),(0,4),(1,5),(2,6),(3,7)]

for row in range(rows):
    for col in range(cols):
        if constallation[row, col] == 0:
            continue  # skip empty cells
        translation = np.array([width * col, length * row, 0])
        # translation defines the staring point to draw the magnet.
        for edge in edges:
            start, end = edge
            start_pt = local_corners_hom[start][:3] + translation
            end_pt = local_corners_hom[end][:3] + translation
            ax.plot([start_pt[0]*1000, end_pt[0]*1000],
                    [start_pt[1]*1000, end_pt[1]*1000],
                    [start_pt[2]*1000, end_pt[2]*1000],
                    color='black', linewidth=1.5)

ax.set_xlabel('x (mm)')
ax.set_ylabel('y (mm)')
ax.set_zlabel('z (mm)')
ax.set_title('B field in 3D (Bx, By, Bz arrows)')


plt.tight_layout()
plt.show()




