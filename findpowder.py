"""
Nov 11, 2024 Developed by Rui Liu @ANL
rui.liu@anl.gov
This module is for detector callibration at 26ID at APS
"""
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from qbin_pilatus import qbin_pilatus
def findpowder(x, img_exp, plotflag=True):
    # Unpack the parameters from x
    R, Tthdet, Gamdet, Xdet, Ydet = x

    # Call qbin_pilatus with the unpacked parameters
    f1 = qbin_pilatus(img_exp, Tthdet, R, Gamdet, Xdet, Ydet, 0, 0, 0, plotflag)

    # If plotflag is True, display the images
    if plotflag:
        # plot the experiment image
        plt.figure(700)
        plt.imshow(img_exp, vmin = 0, vmax = 20, cmap='gray')
        plt.title("Raw Experimental Image")
        plt.colorbar()
        plt.axis('image')

        # Option 1: Plot the Z-component
        plt.figure(701)
        plt.imshow(f1[:, :, 2], cmap='gray')
        plt.title("Z-Component of Data")
        plt.colorbar()
        plt.axis('image')

        # Option 2: Plot the magnitude
        magnitude = np.linalg.norm(f1, axis=2)
        plt.figure(702)
        plt.imshow(magnitude, cmap='gray')
        plt.title("Magnitude of Vector Field")
        plt.colorbar()
        plt.axis('image')

        plt.show()

    # Return the data or error
    if plotflag:
        return f1
    else:
        # Return an error value or some derived quantity
        return np.linalg.norm(f1, axis=2).mean()  # Example: mean of the magnitude
    
# Example usage
# Load your experimental image (replace with your actual image file path)
# img_exp = np.random.rand(1024, 1062)  # Example placeholder image
image_path = 'Images/475_0_25.tiff'
img_exp = imageio.imread(image_path)  # Replace with your image file path
x = [1.8e5, -0.375e5, -0.30e5, 0.00037 * 1e5, 36.97]  # Example parameters

# Call findpowder
result = findpowder(x, img_exp, plotflag=True)