"""
Nov 11, 2024 Developed by Rui Liu @ANL
rui.liu@anl.gov
This module is for detector callibration at 26ID at APS
"""
import numpy as np
import matplotlib.pyplot as plt
from qbin_pilatus import qbin_pilatus
def findpowder(x, img_exp, plotflag=True):
    """
    Main function that uses qbin_pilatus to analyze the image.

    Parameters:
    - x: List of parameters [R, Tthdet, Gamdet, Xdet, Ydet]
    - img_exp: Input experimental image (2D array)
    - plotflag: Boolean flag to determine if plots should be shown (default: True)

    Returns:
    - f: Either the full output of qbin_pilatus or the error value, depending on plotflag
    """
    # Unpack the parameters from x
    R, Tthdet, Gamdet, Xdet, Ydet = x

    # Call qbin_pilatus with the unpacked parameters
    f1 = qbin_pilatus(img_exp, Tthdet, R, Gamdet, Xdet, Ydet, 0, 0, 0, plotflag)

    # If plotflag is True, display the images
    if plotflag:
        plt.figure(701)
        plt.imshow(f1['powder'], cmap='gray')
        plt.title("Thresholded Experiment")
        plt.axis('image')

        plt.figure(702)
        plt.imshow(f1['theory'], cmap='gray')
        plt.title("Theory")
        plt.axis('image')
        plt.show()

    # Return either the full data or just the error, based on plotflag
    if plotflag:
        return f1
    else:
        return f1['error']
    
# Example usage
# Load your experimental image (replace with your actual image file path)
img_exp = np.random.rand(1024, 1062)  # Example placeholder image
x = [1.8e5, -0.375e5, -0.30e5, 0.00037 * 1e5, 36.97]  # Example parameters

# Call findpowder
result = findpowder(x, img_exp, plotflag=True)