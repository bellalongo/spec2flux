# Necessary imports
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.signal import find_peaks
import astropy.units as u
import numpy as np
import pandas as pd
import sys
from os.path import exists
from datetime import date
from collections import defaultdict
from astropy.table import Table
from flux_calc import *
from emission_lines import *
import os
from scipy.integrate import simps

# MAYBE ADD DEF MAIN!!!

# Pull spectra information
filename = sys.argv[1]
instrument = sys.argv[2]
grating = sys.argv[3]
star_name = sys.argv[4]
date = str(date.today())

# Fetch data
grating = grating.upper()
star_name = star_name.upper()
instrument = instrument.upper()
data = fits.getdata(filename)
w, f , e = data['WAVELENGTH'], data['FLUX'], data['ERROR']
mask = (w > 1160) # change if the spectra starts at a different wavelength
wavelength_data, flux_data, error_data = w[mask], f[mask], e[mask]
fresh_start = False # will delete all existing files for that star

# Load Rest Lam data
data = pd.read_csv("DEM_goodlinelist.csv")
rest_lam_data = pd.DataFrame(data)

# Filenames
doppler_filename = "./doppler/" + star_name + "_doppler.txt"
noise_filename = "./noise/" + star_name + "_noise.txt"
emission_line_filename = "./emission_lines" + star_name + "_lines.txt"

# Group emission lines
grouped_lines = grouping_emission_lines(1160, rest_lam_data)

# Find the average width of the peaks
peak_width, peak_width_pixels, flux_range = peak_width_finder(grating, wavelength_data)

# Delete existing data if fresh start
if fresh_start:
    if any(exists([doppler_filename, noise_filename, emission_line_filename])):
        os.remove(doppler_filename)
        os.remove(noise_filename)
        os.remove(emission_line_filename)

# Check if doppler file exists 
doppler_found = exists(doppler_filename)
if doppler_found:
    doppler_shift = np.loadtxt(doppler_filename)*(u.km/u.s)
else:
    doppler_shift = doppler_shift_calc(grouped_lines, wavelength_data, flux_data, flux_range, peak_width, doppler_filename)

# Check if emission line file exists
emission_line_found = exists(emission_line_filename)
if emission_line_found:
    # Load class data! -> UPDATE ME! and go straight to making the final plot ! (should store flux calc? maybe)
    # emission_line_list = loaded data
    print("wow")
# Calculate emission lines
else:
    for index, line in enumerate(emission_line_list):
        # Define necessary variables
        group_wavelength_data = wavelength_data[line.flux_mask]
        group_flux_data = flux_data[line.flux_mask]

        # Necessary emission line info
        rest_lam = line.wavelength_group[len(line.wavelength_group) - 1] * u.AA
        line.obs_lam = doppler_shift.to(u.AA,  equivalencies=u.doppler_optical(rest_lam))
        w0,w1 = wavelength_edges(group_wavelength_data)
        sumerror = (np.sum(error_data[line.flux_mask]**2 * (w1-w0)**2))**0.5

        # If was used for doppler calculation, then not noise
        if line.doppler_candidate:
            noise_bool_list.append(False)
            
            # Calculate continuum and total flux
            continuum = [min(line.fitted_model(group_wavelength_data)) for _ in range(len(group_wavelength_data))]
            total_sumflux = np.sum((line.fitted_model(group_wavelength_data))*(w1-w0)) 
        else:
            # Determine if noise plot
            sns.set_style("darkgrid")
            sns.set_theme(rc={'axes.facecolor':'#F8F5F2'})
            fig = plt.figure(figsize=(14,7))
            ax = fig.add_subplot()
            plt.title(f"Flux vs Wavelength for {line.ion}")
            fig.suptitle("Click 'y' if is noise, 'n' if not", fontweight='bold')
            plt.xlabel('Wavelength (Å)', fontsize =12)
            plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=12)
            cid = fig.canvas.mpl_connect('key_press_event', lambda event: on_key(event, 'Noise Detection'))

            # Check if a fit was successful
            if line.fitted_model:
                voigt_fit, = plt.plot(group_wavelength_data, line.fitted_model(group_wavelength_data), color = "#111D4A") # Change color!
                continuum = [min(line.fitted_model(group_wavelength_data)) for _ in range(len(group_wavelength_data))]
                total_sumflux = np.sum((line.fitted_model(group_wavelength_data))*(w1-w0)) 
            else:
                continuum = split_create_trendline(group_wavelength_data, group_flux_data, peak_width_pixels)
                total_sumflux = np.sum(group_flux_data*(w1-w0))

            # Plot info
            plt.plot(group_wavelength_data, group_flux_data, linewidth=1, color = '#4B3C30')
            continuum_fit, = plt.plot(group_wavelength_data, continuum, color = "#DA667B")

            # Plot all rest wavelengths
            for wavelength in line.wavelength_group:
                rest_lam = plt.axvline(x=wavelength, color = "#71816D", linewidth = 1, linestyle=((0, (5, 5))))
            obs_lam = plt.axvline(x = line.obs_lam.value, color = "#D7816A", linewidth = 1)

            # Plot legend
            if line.fitted_model:
                legend = plt.legend([rest_lam, obs_lam, voigt_fit, continuum_fit], ["Rest Wavelength", "Observed Wavelength", "Voigt Fit", "Continuum"])
            else:
                legend = plt.legend([rest_lam, obs_lam, continuum_fit], ["Rest Wavelength", "Observed Wavelength", "Continuum"])
                
            legend.get_frame().set_facecolor('white')
            plt.show()

        # Check if noise
        if noise_bool_list[index]:
            total_flux = sumerror * (-3)
            sumerror = 0
        else:
            continuum_sumflux = np.sum(continuum*(w1-w0))
            total_flux = total_sumflux - continuum_sumflux
        
        # Update parameters
        line.noise_bool = noise_bool_list[index]    
        line.continuum = continuum
        line.update_flux_error(total_flux, sumerror)

        print(f"Line: {line.obs_lam}, Flux: {line.flux_error[0]}, Error: {line.flux_error[1]}")

# # Create a basic plot
# sns.set_theme()
# plt.title("Flux vs Wavelength for " + star_name)
# plt.xlabel('Wavelength (\AA)')
# plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')
# plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')

# # Plot emission lines
# for line in emission_line_list:
sys.exit()







# Initializing necessary variables
flux = defaultdict(list)
count = 0 
iterations = 0
previous_obs = 0 *u.AA
prev_blended_bool = False
prev_left_bound = 0
emission_lines_list = []


if not noise_found:
    count = 0
    # Determine if the current emission line is noise
    for line in emission_lines_list:
        # Create basic plot
        fig = plt.figure(figsize=(14,7))
        ax = fig.add_subplot()
        fig.suptitle("Click 'y' if should be used for doppler calculation, 'n' if not", fontweight='bold')
        plt.title("Flux vs Wavelength for " + star_name, fontsize=18)
        plt.xlabel('Wavelength (Å)', fontsize=12)
        plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=12)
        trendline_patch = patches.Patch(color='pink', alpha=0.8, label='Continuum')
        rest_patch = patches.Patch(color='lightcoral', alpha=0.8, label='Rest Wavelength')
        obs_patch = patches.Patch(color='darkred', alpha=0.8, label='Observable Wavelength')
        gauss_patch = patches.Patch(color='mediumvioletred', alpha=0.8, label='Gaussian Fit')

        # Plot Gaussian fit
        ax.plot(line.gaussian_x, line.gaussian_y, '-', color='mediumvioletred', linewidth= 2.0)

        # Find Gaussian continuum
        continuum = []
        continuum_array = gaussian_trendline(w[line.flux_mask], line.gaussian_x, line.gaussian_y)

        # Plot emission lines
        ax.plot(w[line.flux_mask], f[line.flux_mask], linewidth = 1.2, alpha = 0.8)
        ax.plot(w[line.flux_mask], continuum_array, color="pink", alpha=0.7)
        plt.axvline(x = line.wavelength, color = 'lightcoral', label = 'Rest wavelength', linewidth= 1.8, ls = '--')
        plt.axvline(x = line.obs_lam.value, color = 'darkred', label = 'Observed wavelength', linewidth= 1.8, ls = '--')
        cid = fig.canvas.mpl_connect('key_press_event', lambda event: on_key(event, 'Noise Detection'))
        plt.legend(handles=[rest_patch, obs_patch, gauss_patch, trendline_patch])
        plt.show()
        
        # Calculate the flux and error
        w0,w1 = wavelength_edges(w[line.flux_mask])
        x_min = np.min(w[line.flux_mask])
        x_max = np.max(w[line.flux_mask])
        total_sumflux = gaussian_integral(amp, mu, sigma, x_min, x_max)
        sumerror = (np.sum(e[line.flux_mask]**2 * (w1-w0)**2))**0.5

        # Calculate the continuum
        for i in range(0, len(continuum_array)):
            continuum.append(continuum_array[i])
        continuum_sumflux = np.sum(continuum*(w1-w0))

        # Check if noise
        total_flux = total_sumflux - continuum_sumflux
        if noise_bool_list[count]:
            # Update emission line's noise bool
            line.noise_bool = noise_bool_list[count]
            # Update flux calculation
            total_flux = sumerror * (-3)
            sumerror = 0

        # Append to flux list
        flux[line.ion].append((line.wavelength, total_flux, sumerror, line.blended_bool))

        count+=1 
    
    # Save to file
    noise_array = np.array(noise_bool_list)
    np.savetxt(noise_filename, noise_array)

    # Printing
    for ion in flux:
        print(f"Ion: {ion} ")
        for data in flux[ion]:
            print(data)

# Plot the emission lines and trendlines
plt.figure(figsize=(14,10))
plt.plot(w[mask], f[mask], color="steelblue")
count = 0

for line in emission_lines_list:
    continuum_array = gaussian_trendline(w[line.flux_mask], line.gaussian_x, line.gaussian_y)

    if line.noise_bool:
        line_color = 'darkgreen'
    else:
        line_color = 'yellowgreen'

    plt.axvline(x=line.obs_lam.value, color= line_color, alpha=0.5)
    trendline, = plt.plot(w[line.flux_mask], continuum_array, color="darkorange", alpha=0.7)

    count+=1

# Create basic plot
plt.title("Flux vs Wavelength for " + star_name)
plt.xlabel('Wavelength (\AA)')
plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')
plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')

# Create plot legend
emission_patch = patches.Patch(color='yellowgreen', alpha=0.7, label='Emission Line')
noise_patch = patches.Patch(color='darkgreen', alpha=0.5, label='Noise')
trendline_patch = patches.Patch(color='darkorange', alpha=0.5, label='Flux Trendline')
plt.legend(handles=[emission_patch, noise_patch, trendline_patch])
plt.show()

# Exit if flux has already been calculated <- NOTE: to adjust flux calculations, DELETE all flux calculations and noise bool
if noise_found:
    sys.exit('Flux already added to ECSV and FITS file, to recalculate, delete noise file and restart')

# Create a fits file
data_array = []
fits_filename = "./flux/" + star_name.lower() + ".fits"
ecsv_filename = "./flux/" + star_name.lower() + ".ecsv"

for ion in flux:
    for data in flux[ion]:
        data_array.append({"Ion": ion, "Wavelength": data[0], "Flux": data[1], "Error": data[2], "Blended line": data[3]})

t = Table(rows=data_array)
t.write(fits_filename, overwrite=True) 
t.write(ecsv_filename, overwrite = True) # does not have a header

# Update header
with fits.open(fits_filename, mode='update') as hdul:
    hdr = hdul[0].header
    hdr.set('DATE', date, 'date flux was calculated')
    hdr.set('FILENAME', filename, 'name of the fits file used to calculate the flux')
    hdr.set('FILETYPE', "SCI", 'file type of fits file')
    hdr.set('TELESCP', "HST", 'telescope used to measure flux')
    hdr.set('INSTRMNT', instrument, 'active instrument to measure flux')
    hdr.set('GRATING', grating, 'grating used to measure flux')
    hdr.set('TARGNAME', star_name, 'name of star used in measurement')
    hdr.set('DOPPLER', str(doppler_shift.value) + " km/s", 'doppler shift used to measure flux')
    hdr.set('WIDTH', "+/- " + str(peak_width) + " Angstroms", 'peak_width used to measure flux')
    hdr.set('RANGE', "+/- " + str(flux_range) + " Angstroms", 'flux range used to measure flux')
    hdr.set('WIDTHPXL', peak_width_pixels, 'peak_width in pixels used to measure flux')
    hdr.set('UPRLIMIT', "3*error", 'upper limit used to determine noise')
    hdul.flush() 