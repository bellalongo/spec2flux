# Necessary imports
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import pandas as pd
import sys
from os.path import exists
from datetime import date
from astropy.table import Table
import os
import json

from emission_lines import *
from flux_calc import *


def main():
    try: 
        filename = sys.argv[1]
        instrument = sys.argv[2]
        grating = sys.argv[3]
        star_name = sys.argv[4]
    except Exception as e:
        sys.exit("Format should be: python main.py spectra/fits_file.fits 'instrument' 'grating' 'star name'")
    today_date = str(date.today())

    # Fetch data
    grating = grating.upper()
    star_name = star_name.upper()
    instrument = instrument.upper()
    data = fits.getdata(filename)
    w, f , e = data['WAVELENGTH'], data['FLUX'], data['ERROR']

    # Adjustable values
    mask = (w > 1160) # change if the spectra starts at a different wavelength
    fresh_start = True # will delete all existing files for that star (set to False if want to just see final plot)
    gaussian_smoothing = False # will smooth spectrum by a preset normalizing term

    # Flux data
    wavelength_data, flux_data, error_data = w[mask], f[mask], e[mask]

    # Find the average width of the peaks
    peak_width, peak_width_pixels, flux_range = peak_width_finder(grating, wavelength_data)

    # Smooth data if gaussian_smoothing
    if gaussian_smoothing:
        wavelength_data, flux_data, error_data = smooth_data(wavelength_data, flux_data, error_data, 1) # adjust sigma as needed

    # Load Rest Lam data
    data = pd.read_csv("DEM_goodlinelist.csv")
    rest_lam_data = pd.DataFrame(data)

    # Load airglow
    airglow_data = pd.read_csv("airglow.csv") 
    airglow_df = pd.DataFrame(airglow_data)

    # Group emission lines 
    grouped_lines = grouping_emission_lines(1160, rest_lam_data)

    # Filenames
    doppler_filename = "./doppler/" + star_name + "_doppler.txt"
    emission_line_filename = "./emission_lines/" + star_name + "_lines.json"
    fits_filename = "./flux/" + star_name.lower() + ".fits"
    ecsv_filename = "./flux/" + star_name.lower() + ".ecsv"
    final_plot_filename = "./plots/" + star_name.lower() + "_final_plot.png"

    # Check run type
    if not fresh_start:
        # Load doppler
        doppler_shift = np.loadtxt(doppler_filename)*(u.km/u.s)

        # Read emission_line data from JSON file
        with open(emission_line_filename, "r") as json_file:
            emission_line_data = json.load(json_file)

        # Reconstruct emission_line objects from dictionaries
        emission_line_list.clear() 
        for data in emission_line_data:
            emission_line_obj = emission_line(**data)
            emission_line_list.append(emission_line_obj)
    else:
        if exists(doppler_filename):
            # Remove all files
            os.remove(doppler_filename)
            os.remove(emission_line_filename)
            os.remove(fits_filename)
            os.remove(ecsv_filename)
            os.remove(final_plot_filename)

        # Calculate doppler
        doppler_shift = doppler_shift_calc(grouped_lines, wavelength_data, flux_data, peak_width, doppler_filename)

        emission_line_data = []
        
        # Iterate through eache emission line
        for index, line in enumerate(emission_line_list):
            # Define emission line variables
            group = line.wavelength_group
            flux_mask = (wavelength_data > group[0] - peak_width) & (wavelength_data < group[len(group) - 1] + peak_width)
            group_wavelength_data = wavelength_data[flux_mask]
            group_flux_data = flux_data[flux_mask]
            group_error_data = error_data[flux_mask]
            group_airglow = airglow_df[(airglow_df['Central Wavelength'] >= np.min(group_wavelength_data)) & (airglow_df['Central Wavelength'] <= np.max(group_wavelength_data))]

            # Emission line information
            rest_lam = line.wavelength_group[len(line.wavelength_group) - 1] * u.AA
            line.obs_lam = doppler_shift.to(u.AA,  equivalencies=u.doppler_optical(rest_lam)).value
            w0,w1 = wavelength_edges(group_wavelength_data)
            sumerror = (np.sum(error_data[flux_mask]**2 * (w1-w0)**2))**0.5

            # If was used for doppler calculation, then not noise
            if line.doppler_candidate:
                noise_bool_list.append(False)
                
                # Calculate continuum and total flux from voigt fit
                continuum = [min(line.fitted_model(group_wavelength_data)) for _ in range(len(group_wavelength_data))]
                total_sumflux = np.sum((line.fitted_model(group_wavelength_data))*(w1-w0)) 
            else:
                # Plot legend parameters and strings
                legend_params, legend_strings = [], []

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
                    voigt_fit, = plt.plot(group_wavelength_data, line.fitted_model(group_wavelength_data), color = "#231651") 
                    continuum = [min(line.fitted_model(group_wavelength_data)) for _ in range(len(group_wavelength_data))]
                    total_sumflux = np.sum((line.fitted_model(group_wavelength_data))*(w1-w0)) 
                    legend_params.append(voigt_fit), legend_strings.append("Voigt Fit")
                else:
                    continuum = split_create_trendline(group_wavelength_data, group_flux_data, peak_width_pixels)
                    total_sumflux = np.sum(group_flux_data*(w1-w0))

                # Plot airglow if applicable
                if len(group_airglow) > 0:
                    for airglow in group_airglow['Central Wavelength']:
                        airglow_lam = plt.axvline(x = airglow, color = '#4464AD', linewidth = 1)
                        legend_params.insert(0, airglow_lam), legend_strings.insert(0, "Airglow")

                # Plot emission lines
                plt.plot(group_wavelength_data, group_flux_data, linewidth=1, color = '#4B3C30')
                continuum_fit, = plt.plot(group_wavelength_data, continuum, color = "#DA667B")
                legend_params.insert(0, continuum_fit), legend_strings.insert(0, "Continuum")
                for wavelength in line.wavelength_group:
                    rest_lam = plt.axvline(x=wavelength, color = "#71816D", linewidth = 1, linestyle=((0, (5, 5))))
                obs_lam = plt.axvline(x = line.obs_lam, color = "#D7816A", linewidth = 1)
                legend_params.insert(0, obs_lam), legend_strings.insert(0, "Observed Wavelength")
                legend_params.insert(0, rest_lam), legend_strings.insert(0, "Rest Wavelength")
                legend = plt.legend(legend_params, legend_strings)
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
            line.continuum = continuum[0] # only need to store 1 value
            line.update_flux_error(total_flux, sumerror)
            emission_line_data.append(emission_line_to_dict(line))

    # Save emission line list information to json
    with open(emission_line_filename, "w") as json_file:
        json.dump(emission_line_data, json_file, indent=4)

    # Create a basic plot
    sns.set_style("darkgrid")
    sns.set_theme(rc={'axes.facecolor':'#F5F5F5'})
    fig = plt.figure(figsize=(14,7))
    ax = fig.add_subplot()
    plt.title("Flux vs Wavelength for " + star_name)
    plt.xlabel('Wavelength (Å)')
    plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
    plt.plot(wavelength_data, flux_data, color = "#0D3B66")

    # Plot emission lines
    for line in emission_line_list:
        # Create mask and continuum
        group = line.wavelength_group
        flux_mask = (wavelength_data > group[0] - peak_width) & (wavelength_data < group[len(group) - 1] + peak_width) # ADJUST ME IF ADJUST FLUX MASK CALC IN EL.PY!
        continuum = [line.continuum for _ in range(len(wavelength_data[flux_mask]))]

        # Check if noise
        if line.noise_bool:
            # Plot rest wavelengths
            for curr_rest in line.wavelength_group:
                noisy_rest_lam = plt.axvline(x=curr_rest, color = '#6F96AE', linewidth = 1.5, linestyle=((0, (5, 5))))
        else:
            # Plot rest and observed wavelengths
            for curr_rest in line.wavelength_group:
                rest_lam = plt.axvline(x=curr_rest, color = '#F3BAD3', linewidth = 1.5, linestyle=((0, (5, 5))))
            obs_lam = plt.axvline(x = line.obs_lam, color = '#DE639A', linewidth = 1.5)

        trendline, = plt.plot(wavelength_data[flux_mask], continuum, color="#EB6424")

    # Plot legend
    legend = plt.legend([noisy_rest_lam, rest_lam, obs_lam, trendline], 
                        ["Noise Wavelength", "Rest Wavelength", "Observed Wavelength", "Continuum"])
    plt.savefig(final_plot_filename)
    plt.show()

    # Check if calculations have already been stored
    if exists(ecsv_filename):
        sys.exit(f"Flux calculations already added to a ecsv and fits file for {star_name}")

    # Add emission line data to an array
    data_array = []
    for line in emission_line_list:
        data_array.append({"Ion": line.ion, "Rest Wavelength" : line.wavelength_group[len(line.wavelength_group) - 1], "Flux" : line.flux_error[0], "Error" : line.flux_error[1], "Blended Line" : line.blended_bool})
        
    # Create file header (EDIT INFORMATION AS NEEDED)
    data_header = [("DATE", today_date, "date flux was calculated"),
                    ("FILENAME", filename, "name of the fits file used to for flux calc"),
                    ("FILETYPE", "SCI", "file type of fits input file"),
                    ("TELESCP", "HST", "telescope used to measure spectrum"),
                    ("INSTRMNT", instrument, "active instrument to measure spectrum"),
                    ("GRATING", grating, "grating used to measure spectrum"),
                    ("TARGNAME", star_name, "name of star used in measurement"),
                    ("DOPPLER", str(doppler_shift.value) + " km/s", "doppler shift used to measure flux"),
                    ("WIDTH", "+/- " + str(peak_width) + " Angstroms", "average peak width of the emissoin lines"),
                    ("RANGE", "+/- " + str(flux_range) + " Angstroms", "flux range used to isolate emission line"), 
                    ("WIDTHPXL", str(peak_width_pixels), "average emission line peak width in pixels"),
                    ("UPRLIMIT", "3*error", "upper limit used for noise")]

    # Create ecsv array to add to ecsv file
    ecsv_table = Table(rows=data_array)
    for header_row in data_header:
        ecsv_table.meta[header_row[0]] = header_row[1]
    ecsv_table.write(ecsv_filename, overwrite = True, format = 'ascii.ecsv')

    # Create a fits table to store the data into the fits file
    dtype = [('Ion', 'S10'), ('Rest_Wavelength', float), ('Flux', float), ('Error', float), ('Blended_Line', bool)]
    data = np.array([(entry['Ion'], entry['Rest Wavelength'], entry['Flux'], entry['Error'], entry['Blended Line']) 
                    for entry in data_array], dtype=dtype)
    hdu = fits.BinTableHDU(data)
    hdu.writeto(fits_filename, overwrite=True)

    # Fits file header
    with fits.open(fits_filename, mode='update') as hdul:
        hdr = hdul[0].header
        for fits_header in data_header:
            hdr.set(fits_header[0], fits_header[1], comment = fits_header[2])
        hdul.flush() 


if __name__ == '__main__':
    main()
    print('All done!')