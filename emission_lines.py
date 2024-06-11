import matplotlib.pyplot as plt
from flux_calc import *
import sys
import math
import astropy.units as u
from astropy.modeling.models import Voigt1D
from astropy.modeling import fitting
import seaborn as sns
from astropy.modeling.fitting import NonFiniteValueError


noise_bool_list = []
doppler_bool_list = []
emission_line_list = []

class emission_line:
    def __init__(self, wavelength_group, ion, obs_lam, noise_bool, blended_bool, doppler_candidate, fitted_model, continuum, flux_error):
        self.wavelength_group = wavelength_group
        self.ion = ion
        self.obs_lam = obs_lam
        self.noise_bool = noise_bool
        self.blended_bool = blended_bool
        self.doppler_candidate = doppler_candidate
        self.fitted_model = fitted_model
        self.continuum = continuum
        self.flux_error = flux_error 

    def update_fitted_model(self, new_fitted_model):
        self.fitted_model = new_fitted_model

    def update_doppler_candidate(self, is_doppler_candidate):
        self.doppler_candidate = is_doppler_candidate
    
    def update_flux_error(self, new_flux, new_error):
        self.flux_error = [new_flux, new_error]


"""
    WILL ONLY SAVE IMPORTANT DATA (no fitted_model) to reduce memory usage
"""
def emission_line_to_dict(emission_line_obj):
    return {
        "wavelength_group": emission_line_obj.wavelength_group,
        "ion": emission_line_obj.ion,
        "obs_lam": float(emission_line_obj.obs_lam),
        "noise_bool": emission_line_obj.noise_bool,
        "blended_bool": emission_line_obj.blended_bool,
        "doppler_candidate": None,
        "fitted_model": None,
        "continuum": float(emission_line_obj.continuum),
        "flux_error": emission_line_obj.flux_error
    }


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
        peak_width = 5.0
    else:
        peak_width = 0.5
        
    flux_range = 2*peak_width

    # Flux range in pixels calculation
    angstroms_to_pixels = wavelength_data[1] - wavelength_data[0] # NOTE! have to recalculate this number every time
    peak_width_pixels = math.floor(peak_width/angstroms_to_pixels)

    return peak_width, peak_width_pixels, flux_range


"""
    Groups blended emission lines together based off of ion, and a pre-determined tolerance
    (Note: this function assumes the DEM data is strictly increasing)
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
        ion = row["Ion"]
        wavelength = float(row["Wavelength"])

        # Check if wavelength is within wavelength bounds
        if wavelength < min_wavelength: 
            continue

        # Check if ion does not exist in the dictionary
        if ion not in ion_groups:
            ion_groups[ion] = [[wavelength]]
        # Ion is already in the dictionary
        else:
            close_group_found = False
            # Iterate through each group for that ion
            for group in ion_groups[ion]:
                # If the largest value in the group - wavelength is less than the tolerance, add to that group
                if abs(max(group) - wavelength) <= tolerance:
                    group.append(wavelength)
                    close_group_found = True
                    break
            
            # If no close group was found, add to that ion as a seperate group
            if not close_group_found:
                ion_groups[ion].append([wavelength])

    return ion_groups


"""
    Calculates the doppler shift based off of peaks and high liklihood rest lam lines
    Name:       doppler_shift_calc()
    Parameters: 
                grouped_lines: grouped emission lines by rest wavelength
                w: current star's masked wavelength data
                f: current star's masked flux data
                peak_width: average emission line width
                doppler_filename: filename where the doppler shift is stores
    Returns:
                doppler_shift: doppler shift of the spectra
"""
def doppler_shift_calc(grouped_lines, w, f, peak_width, doppler_filename):
    rest_candidates, obs_candidates, emission_line_objs = [], [], []

    # Iterate through groups
    for ion in grouped_lines:
        for group in grouped_lines[ion]:
            voigt_profiles = []
            group_mask = (w > group[0] - peak_width) & (w < group[len(group) - 1] + peak_width) 

            # Check if data is valid
            if not any(f[group_mask]):
                continue
            # Iterate through each wavelength in the group
            for wavelength in group:
                # Initial param guesses
                wavelength_mask = (w > wavelength - peak_width/2) & (w < wavelength + peak_width/2) 
                init_x0 = wavelength
                init_amp = np.max(f[wavelength_mask]) 
                init_fwhm_g = peak_width/5
                init_fwhm_l = peak_width/5

                # Voigt profile
                voigt_profile = Voigt1D(x_0 = init_x0, amplitude_L = init_amp, fwhm_L = init_fwhm_l, fwhm_G = init_fwhm_g)
                voigt_profiles.append(voigt_profile)
            
            # Combine voigt profiles
            composite_model = voigt_profiles[0]
            for voigt_profile in voigt_profiles[1:]:
                composite_model += voigt_profile

            # Update emission line list
            group_emission_line_obj = emission_line(
                wavelength_group = group, 
                ion = ion, 
                obs_lam = None, 
                noise_bool = None, 
                blended_bool = True if len(group) > 1 else False, 
                doppler_candidate = None, 
                fitted_model= None,
                continuum = None,
                flux_error = None)
            emission_line_list.append(group_emission_line_obj)

            # Try to fit the model
            try: 
                fitter = fitting.LevMarLSQFitter()
                fitted_model = fitter(composite_model, w[group_mask], f[group_mask])
                group_emission_line_obj.fitted_model = fitted_model
                emission_line_objs.append(group_emission_line_obj)
            except (RuntimeError, TypeError, NonFiniteValueError):
                continue

            # Plot basics
            sns.set_theme()
            fig = plt.figure(figsize=(14,7), facecolor="white")
            ax = fig.add_subplot()
            plt.title(f"Flux vs Wavelength for {ion}")
            fig.suptitle("Click 'y' if should be used for doppler calculation, 'n' if not", fontweight='bold')
            plt.xlabel('Wavelength (Å)', fontsize =12)
            plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=12)
            cid = fig.canvas.mpl_connect('key_press_event', lambda event: on_key(event, 'Doppler Calculation'))

            # Plotting rest and obs emission lines
            group_rest_candidates, group_obs_candidates = [], []
            for i, wavelength in enumerate(group):
                rest_lam = plt.axvline(x=wavelength, color = "#F96E46", linestyle=((0, (5, 5))), linewidth=1)
                group_rest_candidates.append(wavelength)

                if len(group) == 1:
                    obs_lam = plt.axvline(x= fitted_model.x_0.value, color = "#8F1423", linewidth = 1)
                    group_obs_candidates.append(fitted_model.x_0.value)
                    break
                obs_lam = plt.axvline(x= fitted_model[i].x_0.value, color = "#8F1423", linewidth = 1)
                group_obs_candidates.append(fitted_model[i].x_0.value)
            
            plt.plot(w[group_mask], f[group_mask], linewidth=1)
            voigt_fit, = plt.plot(w[group_mask], fitted_model(w[group_mask]), color = "#111D4A")     
            legend = plt.legend([rest_lam, obs_lam, voigt_fit], ["Rest Wavelength", "Observed Wavelength", "Voigt Fit"])
            legend.get_frame().set_facecolor('white')
            plt.show()

            # Update doppler shift calc arrays
            rest_candidates.append(group_rest_candidates)
            obs_candidates.append(group_obs_candidates)

    assert len(doppler_bool_list) > 0, "You didn't click 'y' and 'n' to choose did you?"
    assert len(doppler_bool_list) == len(rest_candidates), "Did you click out of the figure instead of clicking 'y' or 'n'?" 

    # Calculate doppler shift
    dv = []
    # Iterate through each doppler candidate
    for i, boolean in enumerate(doppler_bool_list):
        group_doppler = []
        
        # If is a doppler candidate, set as one in dict
        if boolean: 
            emission_line_objs[i].update_doppler_candidate(True)

            # Calculate doppler for each emission line in that group
            for j, rest_wavelength in enumerate(rest_candidates[i]):
                u_rest_lam = rest_wavelength * u.AA
                u_obs_lam = obs_candidates[i][j] * u.AA
                group_doppler.append(u_obs_lam.to(u.km/u.s,  equivalencies=u.doppler_optical(u_rest_lam)))
            dv.append(sum(group_doppler)/ len(group_doppler))

    doppler_shift = sum(dv)/len(dv)

    # Store value
    with open(doppler_filename, 'a') as f:
        f.write(f"{doppler_shift.value:.3f}\n")
    
    return doppler_shift


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
        print("Invalid key input, select 'y' or 'n'")

    if purpose == "Noise Detection":
        noise_bool_list.append(event.key == 'y')
        plt.close()

    elif purpose == "Doppler Calculation":
        doppler_bool_list.append(event.key == 'y')
        plt.close()

