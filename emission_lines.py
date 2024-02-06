import matplotlib.pyplot as plt
from flux_calc import *
import sys
import math
import matplotlib.patches as patches
from scipy.optimize import curve_fit
from scipy.stats import norm


noise_bool_list = []
doppler_bool_list = []

class emission_line:
    def __init__(self, wavelength, ion, obs_lam, flux_mask, noise_bool, blended_bool):
        self.wavelength = wavelength
        self.ion = ion
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
    if 'L' in grating:
        peak_width = 3.5
    else:
        peak_width = 0.35
        
    flux_range = 2*peak_width

    # Flux range in pixels calculation
    angstroms_to_pixels = wavelength_data[1] - wavelength_data[0] # NOTE! have to recalculate this number every time
    peak_width_pixels = math.floor(peak_width/angstroms_to_pixels)

    return peak_width, peak_width_pixels, flux_range


"""
    Calculates the doppler shift based off of peaks and high liklihood rest lam lines
    Name:       doppler_shift_calc()
    Parameters: 
                rest_lam_data: rest lam data
                wavelength_data: masked wavelength data from the spectra
                flux_range: range to measure the flux of each peak
    Returns:
                doppler_shift: doppler shift of the spectra
"""
def doppler_shift_calc(rest_lam_data, wavelength_data, flux_data, flux_range, star_name):
    # Get rows with high liklihood and best doppler candidates -> ONE PEAK
    """
        instead -> do a guassian fit and then get the peaks of it?
        INSTEAD of find peaks!
    """
    high_liklihood_df = rest_lam_data[rest_lam_data["Likelihood to measure"] == "High"]
    doppler_candidates = []
    obs_candidates = []
    doppler = []

    for index, row in  high_liklihood_df.iterrows():
        # Initalization
        wavelength = row[1]
        flux_range_mask = (wavelength_data > (wavelength - flux_range)) & (wavelength_data < (wavelength + flux_range))

        # See if a gaussian fit can be made
        try:
            # Set up inital parameter guesses 
            init_amp = np.max(flux_data[flux_range_mask])
            init_params = [init_amp, np.mean(wavelength_data[flux_range_mask]), np.std(wavelength_data[flux_range_mask])]

            # Gaussian fit
            popt, _ = curve_fit(gaussian, wavelength_data[flux_range_mask], flux_data[flux_range_mask], p0=init_params)
            amp, mu, sigma = popt
        except RuntimeError:
            continue

        # Create basic plot
        fig = plt.figure(figsize=(14,7))
        ax = fig.add_subplot()
        fig.suptitle('Click "y" if should be used for doppler calc, "n" if not (click "y" for both if doublets)')
        plt.title("Flux vs Wavelength for " + star_name)
        plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ AA$^{-1}$)')
        plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ AA$^{-1}$)')
        # change colors ! -> make into blues instead of reds or purples! 
        rest_patch = patches.Patch(color='darkred', alpha=0.5, label='Rest Wavelength')
        obs_patch = patches.Patch(color = 'lightcoral', alpha = 0.5, label = 'Observed Wavelength Candidates')

        # Plot emission lines
        print(f"Wavelength: {wavelength}")
        ax.plot(wavelength_data[flux_range_mask], flux_data[flux_range_mask], color="steelblue")
        plt.axvline(x = wavelength, color = 'lightcoral', label = 'Emission Line Candidate')
        cid = fig.canvas.mpl_connect('key_press_event', lambda event: on_key(event, 'Doppler Calculation'))

        # Plot the fitted Gaussian profile
        x = np.linspace(np.min(wavelength_data[flux_range_mask]), np.max(wavelength_data[flux_range_mask]), 1000)
        y = gaussian(x, amp, mu, sigma)
        ax.plot(x, y, '-', color='red')

        plt.legend(handles=[rest_patch, obs_patch])
        plt.show()

    sys.exit()



"""
    Gaussian fit function
    Name:       gaussian()
    Parameters: 
                x: data that will be fit
                amp: amplitude of the gaussian fit
                mu: mean of the gaussian fit
                sigma: standard deviation of gaussian fit
    Returns:
                None
"""
def gaussian(x, amp, mu, sigma):
    return amp * norm.pdf(x, mu, sigma)


"""
    Event function that determines if a key was clicked
    Name:       on_key()
    Parameters: 
                event: key press event
    Returns:
                None
"""
def on_key(event, purpose):
    valid_keys = {'y', 'n'}

    if purpose not in ["Noise Detection", "Doppler Calculation"]:
        sys.exit("Invalid purpose, select 'Noise Detection' or 'Doppler Calculation'")

    if event.key not in valid_keys:
        sys.exit("Invalid key input, select 'y' or 'n'")

    if purpose == "Noise Detection":
        noise_bool_list.append(event.key == 'y')
        plt.close()

    elif purpose == "Doppler Calculation":
        doppler_bool_list.append(event.key == 'y')
        plt.close()

