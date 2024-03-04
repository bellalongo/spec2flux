import matplotlib.pyplot as plt
from flux_calc import *
import sys
import math
import matplotlib.patches as patches
from scipy.optimize import curve_fit
from scipy.stats import norm
import astropy.units as u
from scipy.integrate import quad
from collections import defaultdict
from astropy.modeling.models import Voigt1D
from astropy.modeling import models, fitting


noise_bool_list = []
doppler_bool_list = []

class emission_line:
    def __init__(self, wavelength, ion, obs_lam, flux_mask, noise_bool, blended_bool, gaussian_x, gaussian_y):
        self.wavelength = wavelength
        self.ion = ion
        self.obs_lam = obs_lam
        self.flux_mask = flux_mask
        self.noise_bool = noise_bool
        self.blended_bool = blended_bool
        self.gaussian_x = gaussian_x
        self.gaussian_y = gaussian_y


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
    Groups blended emission lines together based off of ion, and a pre-determined tolerance
    (Note: this function assumes the DEM data is increasing)
    Name:       grouping_emission_lines()
    Parameters: 
                rest_lam_data: dataframe of emission lines
    Returns:
                ion_groups: dictionary with ion name as the key, and  blended groups for that as the value
"""
def grouping_emission_lines(min_wavelength, rest_lam_data):
    # Initialize variables
    tolerance = 10.
    ion_groups = {}
    close_group_found = False

    # Loop through emission lines
    for _, row in rest_lam_data.iterrows():
        # Extract ion name and wavelength
        ion = row["Ion"]
        wavelength = float(row["Wavelength"])

        if wavelength < min_wavelength: 
            continue

        # Check if ion already exists in the dictionary
        if ion not in ion_groups:
            ion_groups[ion] = [[wavelength]]
        else:
            # Reset
            close_group_found = False
            for group in ion_groups[ion]:
                # Check if the largest value in the group - wavelength is less than the tolerance
                if abs(max(group) - wavelength) <= tolerance:
                    group.append(wavelength)
                    close_group_found = True
                    break
            
            # If no close group was found
            if not close_group_found:
                ion_groups[ion].append([wavelength])

    for ion in ion_groups:
        print(f"Ion: {ion} {ion_groups[ion]}")

    return ion_groups


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
def doppler_shift_calc(grouped_lines, w, f, flux_range, peak_width, star_name, doppler_filename):
    # Iterate through groups
    for ion in grouped_lines:
        for group in grouped_lines[ion]:
            voigt_profiles = []
            group_mask = (w > group[0] - peak_width) & (w < group[len(group) - 1] + peak_width)
            for wavelength in group:
                # Intialize parameters
                wavelength_mask = (w > wavelength - peak_width/2) & (w < wavelength + peak_width/2)

                init_x0 = wavelength
                init_amp = np.max(f[wavelength_mask])
                init_fwhm_g = 0.1
                init_fwhm_l = 0.1

                # Voigt distributions
                voigt_profile = Voigt1D(x_0 = init_x0, amplitude_L = init_amp, fwhm_L = init_fwhm_l, fwhm_G = init_fwhm_g)
                voigt_profiles.append(voigt_profile)
            
            # Combine voigt distribitions
            composite_model = voigt_profiles[0]
            for voigt_profile in voigt_profiles[1:]:
                composite_model += voigt_profile

            try: 
                # Fit the model
                fitter = fitting.LevMarLSQFitter()
                fitted_model = fitter(composite_model, w[group_mask], f[group_mask])
            except RuntimeError:
                continue
            
            # Plotting
            fig = plt.figure(figsize=(14,7), facecolor="white")
            ax = fig.add_subplot()
            plt.title(f"Flux vs Wavelength for {ion}")
            plt.xlabel('Wavelength (Å)', fontsize =12)
            plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=12)

            plt.plot(w[group_mask], f[group_mask], )
            plt.plot(w[group_mask], fitted_model(w[group_mask]))     
            plt.vlines(group, ymin=min(f[group_mask]), ymax=max(f[group_mask]), colors='pink', linestyles='solid', alpha=0.7)       
            plt.show()

    sys.exit()

    
    # Implement blended line check using checkinrange!
    for _, row in  high_liklihood_df.iterrows():
        # Initalization
        wavelength = row[1]
        flux_range_mask = (wavelength_data > (wavelength - flux_range)) & (wavelength_data < (wavelength + flux_range))

        # See if a gaussian fit can be made
        try:
            # Inital parameter guesses + Gaussian fit
            init_amp = np.max(flux_data[flux_range_mask])
            init_params = [init_amp, np.mean(wavelength_data[flux_range_mask]), np.std(wavelength_data[flux_range_mask])]
            popt, _ = curve_fit(gaussian, wavelength_data[flux_range_mask], flux_data[flux_range_mask], p0=init_params)
            amp, mu, sigma = popt
        except RuntimeError:
            continue

        # Append values to be used for doppler calc
        rest_candidates.append(wavelength)
        obs_candidates.append(mu)

        # Create basic plot
        fig = plt.figure(figsize=(14, 7), facecolor='white')
        ax = fig.add_subplot()
        fig.suptitle("Click 'y' if should be used for doppler calculation, 'n' if not", fontweight='bold')
        plt.title("Flux vs Wavelength for " + star_name, fontsize=18)
        plt.xlabel('Wavelength (Å)', fontsize=12)
        plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=12)
        rest_patch = patches.Patch(color='plum', alpha=0.8, label='Rest Wavelength')
        obs_patch = patches.Patch(color='mediumorchid', alpha=0.8, label='Observable Wavelength')
        gauss_patch = patches.Patch(color='purple', alpha=0.8, label='Gaussian Fit')

        # Plot emission lines 
        ax.plot(wavelength_data[flux_range_mask], flux_data[flux_range_mask], linewidth = 1.2)
        plt.axvline(x=wavelength, color='plum', label='Emission Line Candidate', linewidth=1.8)
        plt.axvline(x=mu, color='mediumorchid', label='Observable Wavelength Candidate', linewidth=1.8)
        cid = fig.canvas.mpl_connect('key_press_event', lambda event: on_key(event, 'Doppler Calculation'))

        # Plot the fitted Gaussian profile
        x = np.linspace(np.min(wavelength_data[flux_range_mask]), np.max(wavelength_data[flux_range_mask]), 1000)
        y = gaussian(x, amp, mu, sigma)
        ax.plot(x, y, '-', color='purple', linewidth=1.5)
        plt.legend(handles=[rest_patch, obs_patch, gauss_patch])
        plt.show()
    
    # Check
    assert len(rest_candidates) == len(obs_candidates)
    
    # Calculate doppler shift
    for i, wavelength in enumerate(rest_candidates):
        u_rest_lam = wavelength * u.AA
        u_obs_lam = obs_candidates[i] * u.AA
        dv.append(u_obs_lam.to(u.km/u.s,  equivalencies=u.doppler_optical(u_rest_lam)))
    doppler_shift = sum(dv)/len(dv)

    # Store value
    with open(doppler_filename, 'a') as f:
        f.write(f"{doppler_shift.value:.3f}\n")

    return doppler_shift

    
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
    Calculate the integral of the Gaussian function over a specified range
    Name: gaussian_integral()
    Parameters:
        amp: amplitude of the gaussian fit
        mu: mean of the gaussian fit
        sigma: standard deviation of the gaussian fit
        x_min: 'a' value (lower value)
        x_max: 'b' value (upper value)
"""
def gaussian_integral(amp, mu, sigma, x_min, x_max):
    integrand = lambda x: amp * norm.pdf(x, mu, sigma)
    return quad(integrand, x_min, x_max)[0]


"""
    Event function that determines if a key was clicked
    Name:       on_key()
    Parameters: 
                event: key press event
                purpose: either doppler or noise
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

