import numpy as np

from create_magnet_constellation import create_magnet_constellation


def magnetic_nodes_halbach(x, y, z, tau, tau_m, hm, Br):
    """
    Calculate magnetic field using the magnetic nodes method.

    Based on Xu's thesis Equations 3.3-3.5.

    Parameters:
        x, y, z   - Coordinates of the measurement point (m).
                    Origin at centre of bottom surface of Halbach array.
        tau       - Pole pitch (m)
        tau_m     - Z-magnet width (m)
        hm        - Magnet height (m)
        Br        - Remanence (T)

    Returns:
        Bx, By, Bz - Magnetic flux density components (T)
    """
    n_z_magnets = 3
    constellation = create_magnet_constellation(n_z_magnets)
    array_size = constellation.shape[0]

    Bx_total = 0.0
    By_total = 0.0
    Bz_total = 0.0

    magnet_spacing = tau / 2

    for i in range(array_size):
        for j in range(array_size):
            magnet_type = constellation[i, j]

            if magnet_type == 0:
                continue

            # Centre position of this magnet (array centred at origin)
            magnet_center_x = (i - (array_size - 1) / 2) * magnet_spacing
            magnet_center_y = (j - (array_size - 1) / 2) * magnet_spacing
            magnet_center_z = hm / 2

            mag_dir = _map_constellation_to_direction(magnet_type)

            if magnet_type in (1, 2):       # up / down  (Z-magnets)
                wx = tau_m
                wy = tau_m
            elif magnet_type in (3, 4):     # left / right (X-magnets)
                wx = (tau - tau_m) / 2
                wy = tau_m
            elif magnet_type in (5, 6):     # X / O  (diagonal magnets)
                wx = (tau - tau_m) / 2
                wy = (tau - tau_m) / 2
            else:
                continue

            Bx_mag, By_mag, Bz_mag = _calculate_magnet_field(
                x, y, z,
                magnet_center_x, magnet_center_y, magnet_center_z,
                wx, wy, hm, Br, mag_dir
            )

            Bx_total += Bx_mag
            By_total += By_mag
            Bz_total += Bz_mag

    return Bx_total, By_total, Bz_total


def _map_constellation_to_direction(code):
    """Map constellation codes to magnetisation direction strings."""
    return {
        1: 'Z+',   # up
        2: 'Z-',   # down
        3: 'X-',   # left
        4: 'X+',   # right
        5: 'XY+',  # X position – diagonal (+X+Y)
        6: 'XY-',  # O position – diagonal (+X-Y)
    }.get(code, 'Z+')


def _calculate_magnet_field(xm, ym, zm, cx, cy, cz, wx, wy, hz, Br, mag_dir):
    """
    Calculate the magnetic field from a single cuboidal magnet.

    (xm, ym, zm)  - measurement point
    (cx, cy, cz)  - magnet centre
    (wx, wy, hz)  - magnet half-widths in x, y and height in z
    Br            - remanence
    mag_dir       - magnetisation direction string
    """
    # 8 corner positions (rows = corners 0..7, cols = x,y,z)
    corners = np.array([
        [cx - wx/2, cy - wy/2, cz - hz/2],  # 0: bottom-front-left
        [cx + wx/2, cy - wy/2, cz - hz/2],  # 1: bottom-front-right
        [cx - wx/2, cy + wy/2, cz - hz/2],  # 2: bottom-back-left
        [cx + wx/2, cy + wy/2, cz - hz/2],  # 3: bottom-back-right
        [cx - wx/2, cy - wy/2, cz + hz/2],  # 4: top-front-left
        [cx + wx/2, cy - wy/2, cz + hz/2],  # 5: top-front-right
        [cx - wx/2, cy + wy/2, cz + hz/2],  # 6: top-back-left
        [cx + wx/2, cy + wy/2, cz + hz/2],  # 7: top-back-right
    ])

    epsilon = np.zeros(8)

    if mag_dir == 'X+':
        epsilon[[0, 2, 4, 6]] = -Br
        epsilon[[1, 3, 5, 7]] = +Br
    elif mag_dir == 'X-':
        epsilon[[0, 2, 4, 6]] = +Br
        epsilon[[1, 3, 5, 7]] = -Br
    elif mag_dir == 'Y+':
        epsilon[[0, 1, 4, 5]] = -Br
        epsilon[[2, 3, 6, 7]] = +Br
    elif mag_dir == 'Y-':
        epsilon[[0, 1, 4, 5]] = +Br
        epsilon[[2, 3, 6, 7]] = -Br
    elif mag_dir == 'Z+':
        epsilon[[0, 1, 2, 3]] = -Br
        epsilon[[4, 5, 6, 7]] = +Br
    elif mag_dir == 'Z-':
        epsilon[[0, 1, 2, 3]] = +Br
        epsilon[[4, 5, 6, 7]] = -Br
    elif mag_dir == 'XY+':
        epsilon[[0, 2, 4, 6]] -= Br / np.sqrt(2)
        epsilon[[1, 3, 5, 7]] += Br / np.sqrt(2)
        epsilon[[0, 1, 4, 5]] -= Br / np.sqrt(2)
        epsilon[[2, 3, 6, 7]] += Br / np.sqrt(2)
    elif mag_dir == 'XY-':
        epsilon[[0, 2, 4, 6]] -= Br / np.sqrt(2)
        epsilon[[1, 3, 5, 7]] += Br / np.sqrt(2)
        epsilon[[0, 1, 4, 5]] += Br / np.sqrt(2)
        epsilon[[2, 3, 6, 7]] -= Br / np.sqrt(2)

    Bx = By = Bz = 0.0

    for k in range(8):
        xk, yk, zk = corners[k]
        xr = xm - xk
        yr = ym - yk
        zr = zm - zk

        r = np.sqrt(xr**2 + yr**2 + zr**2)
        if r < 1e-10:
            continue

        eps_k = epsilon[k]

        # Equations 3.3–3.5
        Bx += (eps_k / (4 * np.pi)) * np.log(-yr + r)
        By += (eps_k / (4 * np.pi)) * np.log(-xr + r)
        Bz += (eps_k / (4 * np.pi)) * np.arctan2(xr * yr, zr * r)

    return Bx, By, Bz




magnetic_nodes_halbach()