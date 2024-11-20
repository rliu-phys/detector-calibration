"""
Nov 11, 2024 Developed by Rui Liu @ANL
rui.liu@anl.gov
This module is for detector callibration at 26ID at APS
"""
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os
from qbin_pilatus import qbin_pilatus

# Function to save results
def save_results(dataout, output_folder="Calibration Results"):
    """
    Save the results to CSV files in the specified output folder.

    Parameters:
    - dataout: Dictionary containing the results (e.g., gamma, two_theta, etc.)
    - output_folder: Folder where results will be saved
    """
    # Create the folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Save each key in the dictionary to a CSV file
    np.savetxt(os.path.join(output_folder, "gamma.csv"), dataout["gamma"], delimiter=",")
    np.savetxt(os.path.join(output_folder, "gamma_mot.csv"), dataout["gamma"], delimiter=",")
    np.savetxt(os.path.join(output_folder, "twotheta.csv"), dataout["two_theta"], delimiter=",")
    np.savetxt(os.path.join(output_folder, "twotheta_mot.csv"), dataout["two_theta"], delimiter=",")

    print(f"Results saved in folder: {output_folder}")

def findpowder(x, img_exp, plotflag=True):
    """
    Function to find powder diffraction calibration.

    Parameters:
    - x: List of parameters [R, Tthdet, Gamdet, Xdet, Ydet]
    - img_exp: Experimental image (2D array)
    - plotflag: Whether to plot the results

    Returns:
    - f: Dictionary containing results if plotflag=True, otherwise the error value
    """
    # Unpack the parameters from x
    Tthdet, R, Gamdet, Xdet, Ydet = x

    # Call qbin_pilatus with the unpacked parameters
    f1 = qbin_pilatus(img_exp, Tthdet, R, Gamdet, Xdet, Ydet, 0, 0, 0, plotflag)

    # If plotflag is True, display the images
    if plotflag:
         # Plot the thresholded experimental data
        plt.figure(701)
        plt.imshow(f1["powder"], cmap="gray")
        plt.title("Thresholded Experiment")
        plt.axis("image")
        plt.colorbar()
        plt.show()

        # Plot the theoretical data
        plt.figure(702)
        plt.imshow(f1["theory"], cmap="gray")
        plt.title("Theory")
        plt.axis("image")
        plt.colorbar()
        plt.show()

        return f1  # Return the full result dictionary

    # If plotflag is False, return only the error
    return f1["error"]

# Example usage
if __name__ == "__main__":
    # Define the image path
    image_path = 'Images/CAP_1.tiff'

    # Load the experimental image
    img_exp = imageio.imread(image_path)

    # Example parameters [R, Tthdet, Gamdet, Xdet, Ydet]
    x = [20, 150000, 17, -18000, 2500]

    # Call the findpowder function with plotting enabled
    result = findpowder(x, img_exp, plotflag=True)
    # Save the results for comparison
    save_results(result)

    # Output the result for debugging
    # print(result)