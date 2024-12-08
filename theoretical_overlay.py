
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio

def calculate_detector_geometry(R, Tthdet, Gamdet, Xdet, Ydet, img_shape, hpxsz=75, vpxsz=75):
    """
    Calculate the detector geometry and pixel positions in real space.
    """
    vNdet, hNdet = img_shape
    vNcen = vNdet // 2
    hNcen = hNdet // 2
    vimaxis = np.arange(1 - vNcen, vNdet - vNcen + 1)
    himaxis = np.arange(1 - hNcen, hNdet - hNcen + 1)

    # Detector center position
    tthd = np.radians(Tthdet)
    gamd = np.radians(Gamdet)
    Detcen = np.array([
        R * np.sin(tthd) * np.cos(gamd),
        R * np.sin(gamd),
        R * np.cos(tthd) * np.cos(gamd)
    ])

    # Normalize detector center vector
    Detunit = Detcen / np.linalg.norm(Detcen)

    # Horizontal and vertical pixel size in micrometers
    pix_x = np.tile(himaxis * hpxsz, (vNdet, 1)) + Xdet
    pix_y = np.tile(vimaxis[:, np.newaxis] * vpxsz, (1, hNdet)) + Ydet

    return pix_x, pix_y, Detunit

def calculate_q_map(pix_x, pix_y, Detunit, R, Ekev):
    """
    Calculate the momentum transfer map (q-magnitude) based on detector geometry.
    """
    kb = 2 * np.pi * Ekev / 12.39842  # Beam momentum in 1/Å

    # Real-space pixel positions
    kfmat = np.zeros((pix_x.shape[0], pix_x.shape[1], 3))
    kfmat[:, :, 0] = pix_x * Detunit[0]
    kfmat[:, :, 1] = pix_y * Detunit[1]
    kfmat[:, :, 2] = pix_x * Detunit[2] + pix_y * Detunit[2] + R

    # Calculate q vector and magnitude
    qmat = np.sqrt(kfmat[:, :, 0]**2 + kfmat[:, :, 1]**2 + kfmat[:, :, 2]**2) * kb
    return qmat

def simulate_theoretical_image(qmat, lat_constants, tolerances):
    """
    Simulate a binary mask for theoretical powder diffraction reflections.
    """
    binary_mask = np.zeros_like(qmat)
    for lat, tol in zip(lat_constants, tolerances):
        for hkl in [[0, 0, 4], [0, 2, 2], [1, 1, 1]]:  # Common reflections
            q_hkl = 2 * np.pi * np.linalg.norm(np.array(hkl) / lat)
            binary_mask += (np.abs(qmat - q_hkl) < tol)
    binary_mask = binary_mask >= 1  # Convert to binary mask
    return binary_mask

def overlay_theoretical_on_experimental(img_exp, binary_mask):
    """
    Overlay theoretical powder reflections on the experimental image.
    """
    impow = img_exp * (1 - binary_mask) + np.max(img_exp) * binary_mask
    return impow

def main(image_path, R, Tthdet, Gamdet, Xdet, Ydet, Ekev=12.0):
    """
    Main function to process the image and overlay theoretical reflections.
    """
    # Load the experimental image
    img_exp = imageio.imread(image_path)
    img_shape = img_exp.shape

    # Detector geometry and pixel calculations
    pix_x, pix_y, Detunit = calculate_detector_geometry(R, Tthdet, Gamdet, Xdet, Ydet, img_shape)

    # Calculate q-map
    qmat = calculate_q_map(pix_x, pix_y, Detunit, R, Ekev)

    # Simulate theoretical powder diffraction image
    lat_constants = [5.4309, 4.0782]  # Silicon and gold lattice constants
    tolerances = [0.007, 0.014]  # Tolerance for Si and Au reflections
    binary_mask = simulate_theoretical_image(qmat, lat_constants, tolerances)

    # Overlay theoretical reflections on experimental image
    impow = overlay_theoretical_on_experimental(img_exp, binary_mask)

    # Plot the experimental image and overlay
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    im0 = axes[0].imshow(img_exp, cmap='viridis', vmax=100)
    axes[0].set_title("Experimental Image")
    axes[0].axis('off')
    fig.colorbar(im0, ax=axes[0], orientation='vertical', label='Intensity')
    
    im1 = axes[1].imshow(binary_mask, cmap='viridis', vmax=100)
    axes[1].set_title("Theoretical Reflection Pattern")
    axes[1].axis('off')
    fig.colorbar(im1, ax=axes[1], orientation='vertical', label='Intensity')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Example usage
    main(
        image_path="IMages/CAP_1.tiff",
        R=1500000,  # Detector distance in micrometers
        Tthdet=20,  # Two-theta angle in degrees
        Gamdet=17,  # Gamma angle in degrees
        Xdet=-18000,  # Detector x-offset in micrometers
        Ydet=2500    # Detector y-offset in micrometers
    )
