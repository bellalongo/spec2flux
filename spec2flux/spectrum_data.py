import astropy.io.fits as fits
from astropy.table import Table
import csv
from datetime import date
import math
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import exists
import pandas as pd
from scipy.ndimage import gaussian_filter
import seaborn as sns


class SpectrumData(object):
    def __init__(self, spectrum_dir, rest_dir, 
                 observation, telescope,
                 instrument, grating, star_name, 
                 min_wavelength, smooth_data):
        
        self.spectrum_dir = spectrum_dir
        self.observation = observation.upper()
        self.telescope = telescope.upper()
        self.instrument = instrument.upper()
        self.grating = grating.upper()
        self.star_name = star_name.upper()
        self.min_wavelength = min_wavelength
        self.todays_date = str(date.today())

        # Fetch data
        data = fits.getdata(spectrum_dir)
        wavelength_data, flux_data, error_data = data['WAVELENGTH'], data['FLUX'], data['ERROR']

        # Mask wavelengths under min_wavelength
        mask = (wavelength_data > min_wavelength)
        self.wavelength_data, self.flux_data, self.error_data = wavelength_data[mask], flux_data[mask], error_data[mask]

        # Get rest wavelengths
        data = pd.read_csv(rest_dir)
        self.rest_lam_data = pd.DataFrame(data)

        # Load airglow
        airglow_data = pd.read_csv("../spec2flux/airglow.csv") # CHECK ME
        self.airglow_df = pd.DataFrame(airglow_data)

        # Get average peak width
        self.line_width, self.line_width_pixels, self.flux_range = self.peak_width_finder()

        # Check if data needs to be smoothed
        if smooth_data:
            self.wavelength_data, self.flux_data, self.error_data = smooth_data(1)

        # Check if spectrum folders exist
        if not exists('doppler'):
            os.makedirs('doppler')
            os.makedirs('emission_lines')
            os.makedirs('flux')
            os.makedirs('plots')
        
        # Spectrum filenames
        self.doppler_dir = f'doppler/{star_name.lower()}_doppler.txt'
        self.emission_lines_dir = f'emission_lines/{star_name.lower()}_lines.json'
        self.fits_dir = f'flux/{star_name.lower()}.fits'
        self.ecsv_dir = f'flux/{star_name.lower()}.ecsv'
        self.csv_dir = f'flux/{star_name.lower()}.csv'
        self.final_plot_dir = f'plots/{star_name.lower()}_final_plot.png' + star_name.lower() + '_final_plot.png'

        # Doppler shift
        self.doppler_shift = None # will be updated in flux_calculator.py


    def peak_width_finder(self):
        """
            Finds the average width of a peak based off of the grating
            Parameters: 
                        None
            Returns:
                        peak_width: average peak width
                        flux_range: range to measure the flux of each peak 
        """
        # Check grating
        if 'L' in self.grating:
            peak_width = 5.0 # NOTE! adjust me as seems fit!
        else:
            peak_width = 0.5
            
        flux_range = 2*peak_width

        # Flux range in pixels calculation
        angstroms_to_pixels = self.wavelength_data[1] - self.wavelength_data[0] 
        peak_width_pixels = math.floor(peak_width/angstroms_to_pixels)

        return peak_width, peak_width_pixels, flux_range
    

    def smooth_data(self, sigma):
        """
            Smooths data using gaussian_filter
            Parameters: 
                        sigma: standard deviation of the gaussian (controls how much neighboring data points contribute to smoothing)
            Returns:
                        smoothed_wavelength: smoothed wavelength data
                        smoothed_flux: smoothed flux data
                        smoothed_error: smoothed error data
        """
        smoothed_wavelength = gaussian_filter(self.wavelength_data, sigma)
        smoothed_flux = gaussian_filter(self.flux_data, sigma)
        smoothed_error = gaussian_filter(self.error_data, sigma)

        return smoothed_wavelength, smoothed_flux, smoothed_error
    

    def save_data(self, emission_lines):
        """
            Adds spectrum data to all file directories in CSV, ECSV, and FITS formats
            Parameters: 
                        emission_lines: EmissionLines instance
            Returns:
                        None

        """
        # Add emission line data to an array
        data_list = []
        for line in emission_lines.line_list:
            data_list.append({"Ion": line.ion, 
                               "Rest Wavelength" : line.group_lam[len(line.group_lam) - 1], 
                               "Flux" : line.flux_error[0], 
                               "Error" : line.flux_error[1], 
                               "Blended Line" : line.blended_bool})
            
        # Create file header (EDIT INFORMATION AS NEEDED)
        data_header = [("DATE", self.todays_date, "date flux was calculated"),
                        ("FILENAME", self.spectrum_dir, "name of the fits file used to for flux calc"),
                        ("FILETYPE", self.observation, "observation type of fits input file"),
                        ("TELESCP", self.telescope, "telescope used to measure spectrum"),
                        ("INSTRMNT", self.instrument, "active instrument to measure spectrum"),
                        ("GRATING", self.grating, "grating used to measure spectrum"),
                        ("TARGNAME", self.star_name, "name of star used in measurement"),
                        ("DOPPLER", str(self.doppler_shift.value) + " km/s", "doppler shift used to measure flux"),
                        ("WIDTH", "+/- " + str(self.line_width) + " Angstroms", "average peak width of the emissoin lines"),
                        ("RANGE", "+/- " + str(self.flux_range) + " Angstroms", "flux range used to isolate emission line"), 
                        ("WIDTHPXL", str(self.line_width_pixels), "average emission line peak width in pixels"),
                        ("UPRLIMIT", "3*error", "upper limit used for noise")]
        
        # Create ecsv array to add to ecsv file
        ecsv_table = Table(rows = data_list)
        for header_row in data_header:
            ecsv_table.meta[header_row[0]] = header_row[1]
        ecsv_table.write(self.ecsv_dir, overwrite = True, format = 'ascii.ecsv')

        # Create CSV column heads and add data to csv file
        data_cols_header = list(data_list[0].keys())
        with open(self.csv_dir, mode='w', newline='') as file:
            writer = csv.writer(file)

            # Write the file header
            writer.writerow(["Header", "Value", "Description"])
            for row in data_header:
                writer.writerow(row)

            # Write an empty row to separate the header from the data
            writer.writerow([])

            # Write the data column headers
            writer.writerow(data_cols_header)

            # Write the data
            dict_writer = csv.DictWriter(file, fieldnames = data_cols_header)
            for data in data_list:
                dict_writer.writerow(data)

        # Create a fits table to store the data into the fits file
        dtype = [('Ion', 'S10'), ('Rest_Wavelength', float), ('Flux', float), ('Error', float), ('Blended_Line', bool)]
        data = np.array([(entry['Ion'], entry['Rest Wavelength'], entry['Flux'], entry['Error'], entry['Blended Line']) 
                        for entry in data_list], dtype=dtype)
        hdu = fits.BinTableHDU(data)
        hdu.writeto(self.fits_dir, overwrite=True)

        # Fits file header
        with fits.open(self.fits_dir, mode='update') as hdul:
            hdr = hdul[0].header
            for fits_header in data_header:
                hdr.set(fits_header[0], fits_header[1], comment = fits_header[2])
            hdul.flush() 


    def final_spectrum_plot(self, emission_lines, flux_calc):
        """
            Plots the final spectrum, with emission lines marked, as well as their respective continuum's and model fitters
            Parameters: 
                        emission_lines: EmissionLines instance
                        flux_calc: FluxCalculator instance
            Returns:
                        None
        """
        # Find max line peak
        max_peak = max(self.flux_data[self.wavelength_data > 1250])
        min_peak = min(self.flux_data)

        # Create a basic plot
        sns.set_style("whitegrid")
        fig = plt.figure(figsize=(14,7))
        ax = fig.add_subplot()
        plt.title("Flux vs Wavelength for " + self.star_name)
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
        plt.plot(self.wavelength_data, self.flux_data, color = "#A7ADC6")
        plt.ylim(min_peak, max_peak * 1.5)

        # Plot emission lines
        for line in emission_lines.line_list:
            # Create mask and continuum
            group = line.group_lam
            flux_mask = (self.wavelength_data > group[0] - self.line_width) & (
                self.wavelength_data < group[len(group) - 1] + self.line_width) # ADJUST ME IF ADJUST FLUX MASK CALC IN EL.PY!
            continuum = [line.continuum for _ in range(len(self.wavelength_data[flux_mask]))]

            # Plot continuum
            trendline, = plt.plot(self.wavelength_data[flux_mask], continuum, color = "#C05746", lw = 2)

            # Plot model fit if there is one
            if line.model_params:
                # Create a model profile and plot
                model_profile = flux_calc.create_model_profile(line.model_params)

                model_fit, = plt.plot(self.wavelength_data[flux_mask], 
                                      model_profile(self.wavelength_data[flux_mask]), color = "#49506F")
                
            # Plot rest wavelengths
            for curr_rest in line.group_lam:
                if line.noise_bool:
                    noisy_rest_lam = plt.axvline(x = curr_rest, color = '#92B257', linewidth = 1.5, linestyle=((0, (5, 5))))
                else:
                    rest_lam = plt.axvline(x = curr_rest, color = '#5D7A3E', linewidth = 1.5, linestyle=((0, (5, 5))))

        # Plot legend
        plt.legend([noisy_rest_lam, rest_lam, model_fit, trendline], 
                            ["Noise Wavelength", "Rest Wavelength", "Model Profile", "Continuum"])
        plt.savefig(self.final_plot_dir)
        plt.show()