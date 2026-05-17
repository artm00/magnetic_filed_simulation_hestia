import numpy as np
# from magnetic_nodes_halbach import magnetic_nodes_halbach
from create_magnet_constellation import create_magnet_constellation


# --- System Parameters (from Xu's thesis Section 3.3.3) ---
tau = 50.8e-3              # Pole pitch (m) = length z magnet plus width x/y magnet
tau_m = 25.4e-3            # Pole width (m) = length/width of z manget
hm = 0.5 * 0.0254          # Magnet thickness: 0.5 inches to metres
Br = 1.32                  # Remanence (T), typical NdFeB

NR_Z_MAGNETS = 6

# create zero pandas array 200x200 where each element is a 3elemented vector (here we will save the B vectors)
# b_field_superimposed =

# get the array of the constallation
constallation = create_magnet_constellation(NR_Z_MAGNETS)
print(constallation)


# i want to do one magnet with its 8 nodes
# dimentions of the magnet:
length = tau_m
width = tau_m
height = hm

sign = -1

# 8 corner nodes of the cuboid; pt_0 is at the origin
pt_0 = np.array([0,      0,     0      , -sign])
pt_1 = np.array([length, 0,     0      , sign])
pt_2 = np.array([0,      width, 0      , sign])
pt_3 = np.array([length, width, 0      , -sign])
pt_4 = np.array([0,      0,     height , sign])
pt_5 = np.array([length, 0,     height , -sign])
pt_6 = np.array([0,      width, height , -sign])
pt_7 = np.array([length, width, height , sign])

node_positions = np.array([pt_0, pt_1, pt_2, pt_3, pt_4, pt_5, pt_6, pt_7])


# TODO Function that calculated the B vector

def compute_b_field(coords, node_positions, epsilon):
    """
    Calculates the 3D magnetic B-field vector at a specific evaluation point.
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
    # print("------------")
    # print(f"coords: {coords}")
    # print(f"x of coords: {coords[0]}")
    # print(f"x of nodes{node_positions[:,0]}")
    # print(R.shape)
    # print(f"x of R coords - node:\n {R[:,0]}")
    # print(f"k: {k}")
    Bx = k * epsilon / (4 * np.pi) * np.log(-R[:,1] + np.sqrt(R[:,0]**2 + R[:,1]**2 + R[:,2]**2))
    By = k * epsilon / (4 * np.pi) * np.log(-R[:,0] + np.sqrt(R[:,0]**2 + R[:,1]**2 + R[:,2]**2))
    Bz = k * epsilon / (4 * np.pi) * np.arctan((R[:,0]*R[:,1]) / (R[:,2]*np.sqrt(R[:,0]**2 + R[:,1]**2 + R[:,2]**2))) ### !!!! maybe z missing, check in the literature again

   
    B_contributions = np.array([Bx, By, Bz])
    B_vector = np.sum(B_contributions, axis=1)
    # print(f"Rz: {R[:,2]}")
    # print(f"Bx: {Bx}")
    # print(B_vector[0])
    # print(f"These values cancel out: {Bx[0]} and {Bx[4]}, pt0 + pt4: {Bx[0] + Bx[4]}")
    # print(f"They should not cancel out as z values are different: {R[0,2]} and {R[4,2]}")

    return B_vector


print(node_positions.shape)

coord_of_interest = np.array([1,1,1])
print(coord_of_interest.shape)
B = compute_b_field(coords=coord_of_interest, node_positions=node_positions, epsilon=Br)
print(B)


# 3D grid — coarse enough to stay readable
nx, ny, nz = 10, 10, 7
x_coords = np.linspace(-2 * length, 3 * length, nx)
y_coords = np.linspace(-2 * width,  3 * width,  ny)
z_coords = np.linspace(-height,     2 * height,  nz)

XXX, YYY, ZZZ = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')

# Compute B at every grid point
B_volume = np.zeros((nx, ny, nz, 3))
for i, x in enumerate(x_coords):
    for j, y in enumerate(y_coords):
        for k, z in enumerate(z_coords):
            coord = np.array([x, y, z])
            B_volume[i, j, k] = compute_b_field(coords=coord, node_positions=node_positions, epsilon=Br)

# PLOTTING
import matplotlib.pyplot as plt

U = B_volume[:, :, :, 0]
V = B_volume[:, :, :, 1]
W = B_volume[:, :, :, 2]

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection='3d')

ax.quiver(XXX * 1000, YYY * 1000, ZZZ * 1000,
          U, V, W,
          length=4, normalize=True, linewidth=0.6, alpha=0.7)

# magnet cuboid outline
lx, ly, lz = length * 1000, width * 1000, height * 1000
corners = np.array([[0,0,0],[lx,0,0],[0,ly,0],[lx,ly,0],
                    [0,0,lz],[lx,0,lz],[0,ly,lz],[lx,ly,lz]])
edges = [(0,1),(0,2),(1,3),(2,3),(4,5),(4,6),(5,7),(6,7),(0,4),(1,5),(2,6),(3,7)]
for a, b in edges:
    ax.plot(*zip(corners[a], corners[b]), color='red', linewidth=1.5)

ax.set_xlabel('x (mm)')
ax.set_ylabel('y (mm)')
ax.set_zlabel('z (mm)')
ax.set_title('B field in 3D (Bx, By, Bz arrows)')
plt.tight_layout()
plt.show()




