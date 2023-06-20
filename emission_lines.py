import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sys
import math

noise_bool_list = []

class emission_line:
    def __init__(self, wavelength, obs_lam, flux_mask, noise_bool, blended_bool):
        self.wavelength = wavelength
        self.obs_lam = obs_lam
        self.flux_mask = flux_mask
        self.noise_bool = noise_bool
        self.blended_bool = blended_bool


"""
    Finds the average width of a peak based off of the grating
    Name:       peak_width_finder()
    Parameters: 
                grating: the grating of the spectra
                wavelength_data: masked wavelength data from the spectra
    Returns:
                peak_width: average peak width
                flux_range: range to measure the flux of each peak 
"""
def peak_width_finder(grating, wavelength_data):
    # Check grating
    if grating == "G140L":
        peak_width = 3.5
    else:
        peak_width = 0.35
        
    flux_range = 2*peak_width

    # Flux range in pixels calculation
    angstroms_to_pixels = wavelength_data[1] - wavelength_data[0] # NOTE! have to recalculate this number every time
    peak_width_pixels = math.floor(peak_width/angstroms_to_pixels)

    return peak_width, peak_width_pixels, flux_range

"""
    Event function that determines if a key was clicked
    Name:       on_key()
    Parameters: 
                event: key press event
    Returns:
                None
"""
def on_key(event):
    if event.key == 'y':
        noise_bool_list.append(True)
        plt.close()
    elif event.key == 'n':
        noise_bool_list.append(False)
        plt.close()
    else:
        sys.exit("Invalid key input")



