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
    
    print('test')
qbin_pilatus([[0]], 0, 0, 0, 0, 0)