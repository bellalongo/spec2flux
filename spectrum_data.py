import astropy.io.fits as fits
from datetime import date
import math
import pandas as pd


class SpectrumData(object):
    def __init__(self, spectrum_dir, rest_dir, instrument, grating, star_name, min_wavelength):
        self.filename = spectrum_dir
        self.instrument = instrument.upper()
        self.grating = grating.upper()
        self.star_name = star_name.upper()
        self.min_wavelength = min_wavelength
        self.todays_date = str(date.today())

        # Fetch data
        data = fits.get_data(spectrum_dir)
        wavelength_data, flux_data, error_data = data['WAVELENGTH'], data['FLUX'], data['ERROR']

        # Mask wavelengths under min_wavelength
        mask = wavelength_data[wavelength_data > min_wavelength]
        self.wavelength_data, self.flux_data, self.error_data = wavelength_data[mask], flux_data[mask], error_data[mask]

        # Get rest wavelengths
        data = pd.read_csv(rest_dir)
        self.rest_lam_data = pd.DataFrame(data)

        # Load airglow
        airglow_data = pd.read_csv("airglow.csv") 
        self.airglow_df = pd.DataFrame(airglow_data)

        # Get average peak width
        self.line_width, self.line_width_pixels, self.flux_range = self.peak_width_finder(grating, wavelength_data)


    def peak_width_finder(self, grating, wavelength_data):
        """
            Finds the average width of a peak based off of the grating
            Parameters: 
                        grating: the grating of the spectra
                        wavelength_data: masked wavelength data from the spectra
            Returns:
                        peak_width: average peak width
                        flux_range: range to measure the flux of each peak 
        """
        # Check grating
        if 'L' in grating:
            peak_width = 5.0 # NOTE! adjust me as seems fit!
        else:
            peak_width = 0.5
            
        flux_range = 2*peak_width

        # Flux range in pixels calculation
        angstroms_to_pixels = wavelength_data[1] - wavelength_data[0] 
        peak_width_pixels = math.floor(peak_width/angstroms_to_pixels)

        return peak_width, peak_width_pixels, flux_range
    

    def smooth_data(wavelength_data, flux_data, error_data, sigma):
        smoothed_wavelength = gaussian_filter(wavelength_data, sigma)
        smoothed_flux = gaussian_filter(flux_data, sigma)
        smoothed_error = gaussian_filter(error_data, sigma)

        return smoothed_wavelength, smoothed_flux, smoothed_error

    # maybe do saving to csv here? 

        

    
