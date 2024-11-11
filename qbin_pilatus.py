"""
Nov 7, 2024 Developed by Rui Liu @ANL
rui.liu@anl.gov
This module is for detector callibration at 26ID at APS
"""
import numpy as np
import matplotlib.pyplot as plt
import imageio
def qbin_pilatus(img_exp, Tthdet, Rdet, Gamdet, Xdet, Ydet, Samth=None, Samchi=None, Samphi=None, plotflag=True):
    """
    Parameters:
    - img_exp: Input experimental image (2D array)
    - Tthdet: Two theta angle to the center of the detector (in degrees)
    - Rdet: Radius from sample to the center of the detector (in micrometers)
    - Gamdet: Out-of-diffraction-plane angle (in degrees)
    - Xdet, Ydet: Detector center coordinates in micrometers
    - Samth, Samchi, Samphi: Sample rotation angles (in degrees)
    - plotflag: Whether to plot images (default: True)
    """
    # Energy and momentum parameters
    EkeV = 12.0 # beam energy in keV
    kb = 2 * np.pi * EkeV / 12.39842 # beam momentum in 1/A

    # pixel size parameters
    hpxsz = 75 # horizontal pixel size in micrometers
    vpxsz = 75 # vertical pixel size in micrometers

    # Sample rotation angles
    if Samth is None: Samth = 24.9456 # Sample rotation about vertical axis (deg, right-hand-rule)
    if Samchi is None: Samchi = 0.0
    if Samphi is None: Samphi = 0.0
    sth = np.radians(Samth)
    sch = np.radians(Samchi)
    sph = np.radians(Samphi)

    # Convert detector angles to radians
    tthd = np.radians(Tthdet)
    gamd = np.radians(Gamdet)

    # Calculate detector center in real space
    Detcen = np.array([
        Rdet * np.sin(tthd) * np.cos(gamd),
        Rdet * np.sin(gamd),
        Rdet * np.cos(tthd) * np.cos(gamd)
    ])
    Detunit = Detcen / np.linalg.norm(Detcen)

    # Calculate unit vectors for pixel directions
    jvec = np.cross([0,1,0], Detunit)
    jvec /= np.linalg.norm(jvec)
    ivec = np.cross(jvec, Detunit)

    # Sample rotation matrices
    sROT1 = np.array([
        [np.cos(sth), 0, np.sin(sth)],
        [0, 1, 0],
        [-np.sin(sth), 0, np.cos(sth)],
    ])
    sROT2 = np.array([
        [np.cos(sch), -np.sin(sch), 0],
        [np.sin(sch), np.cos(sch), 0],
        [0, 0, 1]
    ])
    sROT3 = np.array([
        [1, 0, 0],
        [0, np.cos(sph), -np.sin(sph)],
        [0, np.sin(sph), np.cos(sph)]
    ])
    sROT = sROT1 @ sROT2 @ sROT3

    # Create grid of pixel positions
    vNdet, hNdet = img_exp.shape
    vNcen, hNcen = vNdet // 2, hNdet //2
    vimaxis = np.arange(1 - vNcen, vNdet - vNcen + 1)
    himaxis = np.arange(1 - hNcen, hNdet - hNcen + 1)

    # Generate pixel arrays
    pix_x = np.tile(himaxis * hpxsz + Xdet, (vNdet, 1))
    pix_y = np.tile((vimaxis * vpxsz + Ydet).reshape(-1,1), (1, hNdet))

    # Calculate real-space pixel positions
    kfmat = np.zeros((vNdet, hNdet, 4))
    kfmat[:, :, 0] = pix_x * jvec[0] + pix_y * ivec[0] + Detcen[0]
    kfmat[:, :, 1] = pix_x * jvec[1] + pix_y * ivec[1] + Detcen[1]
    kfmat[:, :, 2] = pix_x * jvec[2] + pix_y * ivec[2] + Detcen[2]
    kfmat[:, :, 3] = np.linalg.norm(kfmat[:, :, :3], axis=2)

    # Normalize kfmat to get unit vectors
    kfmat[:, :, :3] /= kfmat[:, :, 3, np.newaxis]

    # Calculate initial wavevector and q vectors
    sqz = sROT @ np.array([0, 0, 1])
    prefac = -2 * (kfmat[:, :, 0] * sqz[0] + kfmat[:, :, 1] * sqz[1] + kfmat[:, :, 2] * sqz[2])
    kimat = kfmat[:, :, :3] + prefac[:, :, np.newaxis] * sqz

    qmat = kb * (kfmat[:, :, :3] - kimat)
    qmat_magnitude = np.linalg.norm(qmat, axis=2)

    # Plot results if required
    if plotflag:
        plt.figure()
        plt.imshow(qmat_magnitude, cmap='gray')
        plt.title("Momentum transfer (A^-1) per pixel map")
        plt.colorbar()
        plt.show()

    return qmat
# Example usage
# img_exp = imageio.imread('path_to_image.tif')  # Replace with your image file path
# qmat = qbin_pilatus(img_exp, 35.2028, 169510, 0, 10500, -10046)