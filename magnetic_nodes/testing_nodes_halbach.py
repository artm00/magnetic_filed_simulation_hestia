import numpy as np
from b_field_in_magnet_coords import compute_b_field
from first_order_harmonic import compute_b_harmonic_first_order
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

NR_Z_MAGNETS = 6

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
def origin_of_magnet_array(col, row, cols, rows, width, length):
    # return np.array([width * (col), length * (row), 0]) # sets the origin of array at the bottom left corner of the array
    return np.array([width * (col - cols/2), length * (row - rows/2), 0]) # sets the origin of array at the center of the array


def compute_b_volume(x_coords, y_coords, z_coords,
                     constallation, Trafo_dict, node_positions, Br,
                     length, width):
    rows, cols = constallation.shape
    B_volume = np.zeros((len(x_coords), len(y_coords), len(z_coords), 3))

    for row in range(rows):
        for col in range(cols):
            if constallation[row, col] == 0:
                continue
            sign = 1
            if constallation[row, col] in (1, 2):
                if constallation[row, col] == 1:
                    sign = -1
                T = Trafo_dict["R_My"].copy()
            elif constallation[row, col] in (3, 4):
                if constallation[row, col] == 3:
                    sign = -1
                T = Trafo_dict["R_Mx"].copy()
            elif constallation[row, col] in (5, 6):
                if constallation[row, col] == 5:
                    sign = -1
                T = Trafo_dict["R_Mz"].copy()

            translation = origin_of_magnet_array(col, row, cols, rows, width, length)
            translation = np.append(translation, 1)  # make it homogeneous
            T[:, 3] = translation
            T_inv = np.linalg.inv(T)
            Rotation = T[:3, :3]

            for i, x in enumerate(x_coords):
                for j, y in enumerate(y_coords):
                    for k, z in enumerate(z_coords):
                        coord = np.array([x, y, z, 1])
                        coord_mag = T_inv @ coord
                        coord_mag = coord_mag[:3]

                        B_mag = np.array(compute_b_field(
                            coords=coord_mag,
                            node_positions=node_positions,
                            epsilon=sign * Br))

                        B_volume[i, j, k] += Rotation @ B_mag

    return B_volume


# constellation first — its size drives the grid extents
constallation = create_magnet_constellation(NR_Z_MAGNETS)
print(constallation)
rows, cols = constallation.shape

# 3D grid covering the full constellation footprint + one-cell buffer on each side
x_coords = np.linspace(-width - width * cols/2,  cols/2 * width  + width,  10)
y_coords = np.linspace(-length - length * rows/2, rows/2 * length + length, 10)
z_coords = np.linspace(-height * 2, 3* height,              7)

XXX, YYY, ZZZ = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
B_volume = compute_b_volume(x_coords, y_coords, z_coords,
                            constallation, Trafo_dict, node_positions, Br,
                            length, width)
                        

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
        translation = origin_of_magnet_array(col, row, cols, rows, width, length)
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
ax.set_box_aspect([
    x_coords[-1] - x_coords[0],
    y_coords[-1] - y_coords[0],
    z_coords[-1] - z_coords[0],
])

plt.tight_layout()
plt.show()


# ------------------------- LINE SAMPLE PLOT -------------------------

def compute_b_along_line(points, constallation, Trafo_dict, node_positions, Br, length, width):
    """Sample B field at each point in `points` (shape N×3, world coords)."""
    rows, cols = constallation.shape
    B_line = np.zeros((len(points), 3))

    for row in range(rows):
        for col in range(cols):
            if constallation[row, col] == 0:
                continue
            sign = 1
            if constallation[row, col] in (1, 2):
                if constallation[row, col] == 1:
                    sign = -1
                T = Trafo_dict["R_My"].copy()
            elif constallation[row, col] in (3, 4):
                if constallation[row, col] == 3:
                    sign = -1
                T = Trafo_dict["R_Mx"].copy()
            elif constallation[row, col] in (5, 6):
                if constallation[row, col] == 5:
                    sign = -1
                T = Trafo_dict["R_Mz"].copy()

            translation = origin_of_magnet_array(col, row, cols, rows, width, length)
            translation = np.append(translation, 1)
            T[:, 3] = translation
            T_inv = np.linalg.inv(T)
            Rotation = T[:3, :3]

            for idx, pt in enumerate(points):
                coord = np.append(pt, 1)
                coord_mag = (T_inv @ coord)[:3]
                B_mag = np.array(compute_b_field(
                    coords=coord_mag,
                    node_positions=node_positions,
                    epsilon=sign * Br))
                B_line[idx] += Rotation @ B_mag

    return B_line


def plot_b_along_line(p_start, p_end, n_points, constallation, Trafo_dict,
                      node_positions, Br, length, width, tau, tau_m, hm):
    t = np.linspace(0, 1, n_points)
    points = np.outer(1 - t, p_start) + np.outer(t, p_end)  # (n_points, 3)
    dist_mm = t * np.linalg.norm(p_end - p_start) * 1000

    B_nodes    = compute_b_along_line(points, constallation, Trafo_dict,
                                      node_positions, Br, length, width)
    B_harmonic = compute_b_harmonic_first_order(points, tau, tau_m, hm, Br)

    labels = ['Bx', 'By', 'Bz']
    _, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for k, ax in enumerate(axes):
        ax.plot(dist_mm, B_nodes[:, k],    'b-',  linewidth=1.5, label='Magnetic nodes')
        ax.plot(dist_mm, B_harmonic[:, k], 'r--', linewidth=1.5, label='Harmonic (1st order)')
        ax.set_ylabel(f'{labels[k]} (T)')
        ax.legend(loc='best')
        ax.grid(True)
    axes[-1].set_xlabel('Distance along line (mm)')
    axes[0].set_title(
        f'B field along line  start=({p_start[0]*1e3:.1f}, {p_start[1]*1e3:.1f}, {p_start[2]*1e3:.1f}) mm  '
        f'end=({p_end[0]*1e3:.1f}, {p_end[1]*1e3:.1f}, {p_end[2]*1e3:.1f}) mm')
    plt.tight_layout()
    plt.show()


p_start = np.array([3.5e-3, 10.5e-3, -25e-3])
p_end   = np.array([192.5e-3, 37.5e-3, -25e-3])
plot_b_along_line(p_start, p_end, 200, constallation, Trafo_dict,
                  node_positions, Br, length, width, tau, tau_m, hm)
