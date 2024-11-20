"""
Nov 12, 2024 Developed by Rui Liu @ANL
rui.liu@anl.gov
This module is for detector callibration at 26ID at APS
"""
import numpy as np
import matplotlib.pyplot as plt

def qbin_pilatus(img_exp, Tthdet, Rdet, Gamdet, Xdet, Ydet, Samth=24.9456, Samchi=0.0, Samphi=0.0, plotflag=True):
    """
    Function to perform detector calibration and momentum transfer calculations.
    
    Parameters:
    - img_exp: Experimental image (2D array)
    - Tthdet, Rdet, Gamdet, Xdet, Ydet: Detector parameters
    - Samth, Samchi, Samphi: Sample rotation angles (defaults provided)
    - plotflag: Whether to plot images
    
    Returns:
    - dataout: Dictionary containing calculated data
    """

    # Constants and parameters
    Ekev = 12.0  # Beam energy in keV
    kb = 2 * np.pi * Ekev / 12.39842  # Beam momentum in 1/Å
    hpxsz = 75  # Horizontal pixel size in micrometers
    vpxsz = 75  # Vertical pixel size
    Binflag = 1  # Set to 1 if binning is required
    Powdflag = 1  # Set to 1 to display powder lines
    qnum = 75  # Number of momentum space bins

    # Load the experimental image
    imdet = img_exp.astype(np.float64)  # Convert image to float for calculations
    vNdet, hNdet = imdet.shape
    vNcen = vNdet // 2
    hNcen = hNdet // 2
    vimaxis = np.arange(1 - vNcen, vNdet - vNcen + 1)
    himaxis = np.arange(1 - hNcen, hNdet - hNcen + 1)

    # Detector parameters in radians
    tthd = np.radians(Tthdet)  # Convert two theta to radians
    gamd = np.radians(Gamdet)  # Convert gamma to radians

    # Calculate detector center in real space (micrometers)
    Detcen = np.array([
        Rdet * np.sin(tthd) * np.cos(gamd),
        Rdet * np.sin(gamd),
        Rdet * np.cos(tthd) * np.cos(gamd)
    ])

    # Normalize to get the unit vector
    Detunit = Detcen / np.linalg.norm(Detcen)

    # Calculate the horizontal and vertical unit vectors for the pixel plane
    jvec = np.cross([0, 1, 0], Detunit)  # Cross product to get horizontal vector
    jvec /= np.linalg.norm(jvec)  # Normalize
    ivec = np.cross(jvec, Detunit)  # Cross product to get vertical vector

    # Generate pixel arrays (in micrometers)
    pix_x = np.tile(himaxis * hpxsz, (vNdet, 1)) + Xdet
    pix_y = np.tile(vimaxis[:, np.newaxis] * vpxsz, (1, hNdet)) + Ydet

    # Calculate real-space pixel positions (in micrometers)
    kfmat = np.zeros((vNdet, hNdet, 4))
    kfmat[:, :, 0] = pix_x * jvec[0] + pix_y * ivec[0] + Detcen[0]
    kfmat[:, :, 1] = pix_x * jvec[1] + pix_y * ivec[1] + Detcen[1]
    kfmat[:, :, 2] = pix_x * jvec[2] + pix_y * ivec[2] + Detcen[2]

    # Calculate the magnitude of the position vector
    kfmat[:, :, 3] = np.sqrt(kfmat[:, :, 0]**2 + kfmat[:, :, 1]**2 + kfmat[:, :, 2]**2)

    # Normalize the vectors
    kfmat[:, :, 0] /= kfmat[:, :, 3]
    kfmat[:, :, 1] /= kfmat[:, :, 3]
    kfmat[:, :, 2] /= kfmat[:, :, 3]

    # Sample rotation angles in radians
    sth = np.radians(Samth)
    sch = np.radians(Samchi)
    sph = np.radians(Samphi)

    # Sample rotation matrices
    sROT1 = np.array([
        [np.cos(sth), 0, np.sin(sth)],
        [0, 1, 0],
        [-np.sin(sth), 0, np.cos(sth)]
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

    # Combined rotation matrix
    sROT = sROT1 @ sROT2 @ sROT3
    sqz = sROT @ np.array([1, 0, 0])  # 00L rod direction

    # Calculate kinitial vectors
    prefac = -2 * (kfmat[:, :, 0] * sqz[0] + kfmat[:, :, 1] * sqz[1] + kfmat[:, :, 2] * sqz[2])
    kimat = np.zeros_like(kfmat)
    kimat[:, :, 0] = kfmat[:, :, 0] + prefac * sqz[0]
    kimat[:, :, 1] = kfmat[:, :, 1] + prefac * sqz[1]
    kimat[:, :, 2] = kfmat[:, :, 2] + prefac * sqz[2]

    # Calculate the momentum transfer vector q
    qmat = np.zeros_like(kimat)
    qmat[:, :, 0] = kb * (kfmat[:, :, 0] - kimat[:, :, 0])
    qmat[:, :, 1] = kb * (kfmat[:, :, 1] - kimat[:, :, 1])
    qmat[:, :, 2] = kb * (kfmat[:, :, 2] - kimat[:, :, 2])
    qmat[:, :, 3] = np.sqrt(qmat[:, :, 0]**2 + qmat[:, :, 1]**2 + qmat[:, :, 2]**2)

    # Binning in magnitude q
    qlist = np.zeros((qnum, 2))
    if Binflag:
        # Extract the magnitude of q
        q_values = qmat[:, :, 3]
        qmin = np.min(q_values[q_values > 0])  # Minimum non-zero q value
        qmax = np.max(q_values)
        for i in range(qnum):
            q1 = qmin + (i * (qmax - qmin) / qnum)
            q2 = qmin + ((i + 1) * (qmax - qmin) / qnum)
            qlist[i, 0] = q1

            # Create masks for binning
            mask = (q_values >= q1) & (q_values < q2)
            if np.sum(mask) > 0:
                qlist[i, 1] = np.sum(imdet[mask]) / np.sum(mask)

        # Plot the binned average intensity
        if plotflag:
            plt.figure()
            plt.plot(qlist[:, 0], np.log(qlist[:, 1] + 1))
            plt.xlabel("Momentum transfer (1/Å)")
            plt.ylabel("Average Intensity (ADU)")
            plt.title("Binned Average Intensity")
            plt.show()

    # Display powder lines if Powdflag is set
    if Powdflag:
        tempm1 = np.zeros((vNdet, hNdet))
        iref, jref = [], []
        iref2, jref2 = [], []
        # Si and Au lattice reflections
        refl = np.array([[0, 0, 4], [0, 2, 2], [1, 1, 1], [1, 1, 3], [1, 3, 3], [3, 3, 3], [0, 4, 4]])
        refl2 = np.array([[0, 0, 4], [0, 0, 2], [0, 2, 2], [1, 1, 1], [1, 1, 3], [1, 3, 3], [3, 3, 3], [0, 4, 4]])
        lat = np.array([5.4309, 5.4309, 5.4309])
        lat2 = np.array([4.0782, 4.0782, 4.0782])

        # Loop through Si reflections
        for reflection in refl:
            q1 = 2 * np.pi * np.linalg.norm(reflection / lat)
            mask = np.abs(qmat[:, :, 3] - q1) < 0.007
            tempm1 += mask
            r, c = np.where(mask)
            if r.size > 0:
                iref.append(r[0])
                jref.append(c[0])
            else:
                iref.append(-1)
                jref.append(-1)

        # Loop through Au reflections
        for reflection in refl2:
            q1 = 2 * np.pi * np.linalg.norm(reflection / lat2)
            mask = np.abs(qmat[:, :, 3] - q1) < 0.014
            tempm1 += mask
            r, c = np.where(mask)
            if r.size > 0:
                iref2.append(r[0])
                jref2.append(c[0])
            else:
                iref2.append(-1)
                jref2.append(-1)

        tempm1 = tempm1 >= 1
        impow = imdet * (1 - tempm1) + np.max(imdet) * tempm1

        # Plot the image with powder lines
        plt.figure()
        plt.imshow(impow, cmap='gray')
        plt.axis('image')
        plt.title("Image with Powder Lines")
        for i, (y, x) in enumerate(zip(iref, jref)):
            if y >= 0:
                plt.text(x, y, f"Si[{refl[i, 0]}{refl[i, 1]}{refl[i, 2]}]", fontsize=8, color='white')
        for i, (y, x) in enumerate(zip(iref2, jref2)):
            if y >= 0:
                plt.text(x, y, f"Au[{refl2[i, 0]}{refl2[i, 1]}{refl2[i, 2]}]", fontsize=8, color='white')
        plt.show()

    # Prepare the final data output
    dataout = {
        "two_theta": np.degrees(np.arccos(kfmat[:, :, 2])),
        "gamma": 90 - np.degrees(np.arccos(kfmat[:, :, 1])),
        "qmat": qmat,
        "kfmat": kfmat,
        "kimat": kimat,
        "qlist": qlist
    }

    return dataout
