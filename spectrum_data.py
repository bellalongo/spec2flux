import astropy.io.fits as fits
from datetime import date
import math
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


class SpectrumData(object):
    def __init__(self, spectrum_dir, rest_dir, 
                 instrument, grating, star_name, 
                 min_wavelength, smooth_data):
        
        self.filename = spectrum_dir
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
        airglow_data = pd.read_csv("airglow.csv") 
        self.airglow_df = pd.DataFrame(airglow_data)

        # Get average peak width
        self.line_width, self.line_width_pixels, self.flux_range = self.peak_width_finder()

        # Check if data needs to be smoothed
        if smooth_data:
            self.wavelength_data, self.flux_data, self.error_data = smooth_data(1)
        
        # Spectrum filenames
        self.doppler_dir = './doppler/' + star_name + '_doppler.txt'
        self.emission_lines_dir = './emission_lines/' + star_name + '_lines.json'
        self.fits_dir = './flux/' + star_name.lower() + '.fits' # adjust filenames
        self.ecsv_dir = './flux/' + star_name.lower() + '.ecsv'
        self.csv_dir = './flux/' + star_name.lower() + '.csv' # make it so they're added to a CSV
        self.final_plot_dir = './plots/' + star_name.lower() + '_final_plot.png'


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

    # maybe do saving to csv here? 

        

    
