import astropy.io.fits as fits
from astropy.table import Table
import csv
from datetime import date
import math
import numpy as np
import os
from os.path import exists
import pandas as pd
from scipy.ndimage import gaussian_filter


# ------------------------------
# Constants
# ------------------------------
PEAK_WIDTH_LOW_RES = 5.0  # Angstroms
PEAK_WIDTH_HIGH_RES = 0.5  # Angstroms

# ------------------------------
# SpectrumData Class
# ------------------------------
class SpectrumData(object):
    """
        This class manages the loading of spectral data from FITS or ECSV files,
        applies necessary preprocessing, and provides methods for saving results.
        It handles wavelength masking, smoothing, and maintains file paths for
        storing various outputs of the analysis.
        Attributes:
            spectrum_dir (str): Path to the spectrum file

            observation (str): Observation type (e.g., 'SCI')
            telescope (str): Telescope used (e.g., 'HST')
            instrument (str): Instrument used (e.g., 'STIS', 'COS')
            grating (str): Grating type
            resolution (str): Resolution level ('HIGH' or 'LOW')
            star_name (str): Name of the target star
            min_wavelength (float): Minimum wavelength for analysis
            max_wavelength (float): Maximum wavelength for analysis
            todays_date (str): Current date in string format
            apply_smoothing (bool): Whether to apply Gaussian smoothing
            fresh_start (bool): Whether to start a new analysis

            line_fit_model (str): Model for fitting emission lines
            cont_fit (str): Continuum fitting method

            wavelength_data (ndarray): Array of wavelength values
            flux_data (ndarray): Array of flux values
            error_data (ndarray): Array of error values
            rest_lam_data (DataFrame): Rest wavelength data for emission lines

            airglow_df (DataFrame): Airglow data

            line_width (float): Width of emission lines in Angstroms
            line_width_pixels (int): Width of emission lines in pixels
            flux_range (float): Range of flux for isolating emission lines

            doppler_shift (Quantity): Calculated Doppler shift
            spectrum_continuum (float): Calculated spectrum continuum level
    """
    def __init__(self, spectrum_config, analysis_config):
        """
            Initializes the SpectrumData object with configuration settings.
            Arguments:
                spectrum_config (dict): Dictionary containing spectrum configuration values
                analysis_config (dict): Dictionary containing analysis configuration values
        """
        self.spectrum_config = spectrum_config
        self.analysis_config = analysis_config

        # Load the data and setup data directories
        self._load_data()
        self._setup_directories()

        # Delete files if fresh start
        if self.fresh_start: self._delete_existing_files()

        # Initialize doppler shift as None
        self.doppler_shift = None

        # Initialize spectrum continuum
        self.spectrum_continuum = None

    # ------------------------------
    # Private Helper Functions
    # ------------------------------
    def _load_data(self):
        """
            Loads spectrum data from specified files and applies initial processing.
            Arguments:
                None
            Returns: 
                None
            Raises:
                ValueError: If the spectrum file format is not supported
        """
        # Extract from spectrum config
        self.spectrum_dir = self.spectrum_config['spectrum_dir']
        self.observation = self.spectrum_config['observation'].upper()
        self.telescope = self.spectrum_config['telescope'].upper()
        self.instrument = self.spectrum_config['instrument'].upper()
        self.grating = self.spectrum_config['grating'].upper()
        self.resolution = self.spectrum_config['resolution'].upper()
        self.star_name = self.spectrum_config['star_name'].upper()
        self.min_wavelength = self.spectrum_config['min_wavelength']
        self.max_wavelength = self.spectrum_config['max_wavelength']
        self.todays_date = str(date.today())

        # Extract from analysis config
        self.apply_smoothing = self.analysis_config['apply_smoothing']
        self.line_fit_model = self.analysis_config['line_fit_model']
        self.cont_fit = self.analysis_config['cont_fit']
        self.fresh_start = self.analysis_config['fresh_start']

        # Get data from files
        data = self._load_spectrum_file()
        wavelength_data, flux_data, error_data = self._extract_data(data)
        
        # Apply wavelength mask
        mask = (wavelength_data > self.min_wavelength) & (wavelength_data < self.max_wavelength)
        self.wavelength_data = wavelength_data[mask]
        self.flux_data = flux_data[mask]
        self.error_data = error_data[mask]

        # Load additional data
        self.rest_lam_data = pd.read_csv(self.spectrum_config['rest_dir'])
        self.airglow_df = pd.read_csv(self.spectrum_config['airglow_dir'])

        # Calculate peak width
        self.line_width, self.line_width_pixels, self.flux_range = self.peak_width_finder()

        # Apply smoothing if needed
        if self.apply_smoothing:
            self.wavelength_data, self.flux_data, self.error_data = self.smooth_data(1)

    def _load_spectrum_file(self):
        """
            Loads the spectrum file based on its extension.
            Arguments:
                None
            Returns:
                object: Data structure containing the spectrum data
            Raises:
                ValueError: If the file format is not supported (.fits or .ecsv)
        """
        # Check if fits file
        if self.spectrum_dir.endswith('.fits'):
            return fits.getdata(self.spectrum_dir)
        # Check if escv file
        elif self.spectrum_dir.endswith('.ecsv'):
            return Table.read(self.spectrum_dir, format='ascii.ecsv')
        raise ValueError("File format not supported. Please provide a .fits or .ecsv file.")

    def _extract_data(self, data):
        """
            Extracts wavelength, flux, and error arrays from the loaded data.
            Arguments:
                data: The loaded spectrum data structure
            Returns:
                tuple: (wavelength_data, flux_data, error_data) arrays
            Raises:
                None
        """
        try:
            return data['WAVELENGTH'], data['FLUX'], data['ERROR']
        except KeyError:
            return data[0][0], data[0][1], data[0][2]

    def _setup_directories(self):
        """
            Creates necessary directories for output files and sets file paths.
            Arguments:
                None
            Returns:
                None
            Raises:
                None
        """
        # Initialize directory roots for data to be saved to
        directories = ['doppler', 'emission_lines', 'flux', 'plots']
        for directory in directories:
            if not exists(directory):
                os.makedirs(directory)

        # Create directory paths
        base_name = self.star_name.lower()
        self.doppler_dir = f'doppler/{base_name}_doppler.txt'
        self.emission_lines_dir = f'emission_lines/{base_name}_lines.json'
        self.fits_dir = f'flux/{base_name}.fits'
        self.ecsv_dir = f'flux/{base_name}.ecsv'
        self.csv_dir = f'flux/{base_name}.csv'
        self.final_plot_dir = f'plots/{base_name}_final_plot.png'

    def _prepare_data_list(self, emission_lines):
        """
            Prepares a list of emission line data for output files.
            Arguments:
                emission_lines (EmissionLines): Object containing emission line data
            Returns:
                list: List of dictionaries with emission line information
            Raises:
                None
        """
        data_list = []

        # Iterate through each emission line
        for ion in emission_lines.line_list:
            # Append emission line data 
            for _, line in emission_lines.line_list[ion].items():
                data_list.append({
                    "Ion": line.ion,
                    "Rest Wavelength": line.group_lam[-1],
                    "Flux": line.flux_error[0],
                    "Error": line.flux_error[1],
                    "Blended Line": line.blended_bool,
                })
        return data_list

    def _prepare_data_header(self):
        """
            Prepares header information for output files.
            Arguments:
                None
            Returns:
                list: List of tuples containing header field names, values, and descriptions
            Raises:
                None
        """
        return [
            ("DATE", self.todays_date, "date flux was calculated"),
            ("FILENAME", self.spectrum_dir, "name of the fits file used for flux calc"),
            ("FILETYPE", self.observation, "observation type of fits input file"),
            ("TELESCP", self.telescope, "telescope used to measure spectrum"),
            ("INSTRMNT", self.instrument, "active instrument to measure spectrum"),
            ("GRATING", self.grating, "grating used to measure spectrum"),
            ("TARGNAME", self.star_name, "name of star used in measurement"),
            ("DOPPLER", f"{self.doppler_shift.value} km/s", "doppler shift used to measure flux"),
            ("WIDTH", f"+/- {self.line_width} Angstroms", "average peak width of the emission lines"),
            ("RANGE", f"+/- {self.flux_range} Angstroms", "flux range used to isolate emission line"),
            ("WIDTHPXL", str(self.line_width_pixels), "average emission line peak width in pixels"),
            ("UPRLIMIT", "3*error", "upper limit used for noise"),
        ]

    def _save_ecsv(self, data_list, data_header):
        """
            Saves emission line data as an ECSV file.
            Arguments:
                data_list (list): List of dictionaries with emission line information
                data_header (list): List of header information tuples
            Returns:
                None
            Raises:
                None
        """
        # Grab ECSV table
        ecsv_table = Table(rows=data_list)

        # Irerate through the heading
        for header_row in data_header:
            ecsv_table.meta[header_row[0]] = header_row[1]

        # Save the table to the ECSV path
        ecsv_table.write(self.ecsv_dir, overwrite=True, format='ascii.ecsv')

    def _save_csv(self, data_list, data_header):
        """
            Saves emission line data as a CSV file.
            Arguments:
                data_list (list): List of dictionaries with emission line information
                data_header (list): List of header information tuples
            Returns:
                None
            Raises:
                None
        """
        # Grab the header
        data_cols_header = list(data_list[0].keys())

        with open(self.csv_dir, mode='w', newline='') as file:
            writer = csv.writer(file)

            # Create header values 
            writer.writerow(["Header", "Value", "Description"])

            # Create header 
            for row in data_header:
                writer.writerow(row)
            
            # Gap in header 
            writer.writerow([])
            writer.writerow(data_cols_header)

            # Grab data and write to the csv file
            dict_writer = csv.DictWriter(file, fieldnames=data_cols_header)
            for data in data_list:
                dict_writer.writerow(data)

    def _save_fits(self, data_list, data_header):
        """
            Saves emission line data as a FITS file.
            Arguments:
                data_list (list): List of dictionaries with emission line information
                data_header (list): List of header information tuples
            Returns:
                None
            Raises:
                None
        """
        # Create column description
        dtype = [
            ('Ion', 'S10'),
            ('Rest_Wavelength', float),
            ('Flux', float),
            ('Error', float),
            ('Blended_Line', bool),
        ]
        
        # Initailize column data 
        data = np.array(
            [
                (entry['Ion'], entry['Rest Wavelength'], entry['Flux'], entry['Error'], entry['Blended Line'])
                for entry in data_list
            ],
            dtype=dtype,
        )

        # Save column data 
        hdu = fits.BinTableHDU(data)
        hdu.writeto(self.fits_dir, overwrite=True)

        # Create header 
        with fits.open(self.fits_dir, mode='update') as hdul:
            hdr = hdul[0].header
            for fits_header in data_header:
                hdr.set(fits_header[0], fits_header[1], comment=fits_header[2])
            hdul.flush()

    def _delete_existing_files(self):
        """
            Deletes existing output files if fresh_start is True.
            Arguments:
                None
            Returns:
                None
            Raises:
                None
        """
        # Save all paths asociated with the current star 
        files_to_delete = [
            self.doppler_dir,
            self.emission_lines_dir,
            self.fits_dir,
            self.ecsv_dir,
            self.csv_dir,
            self.final_plot_dir
        ]

        # Iterate through each path
        for file_path in files_to_delete:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
            except OSError as e:
                print(f"Warning: Could not delete {file_path}: {e}")

    # ------------------------------
    # Public Methods
    # ------------------------------
    def peak_width_finder(self):
        """
            Determines the appropriate peak width based on resolution.
            Arguments:
                None
            Returns:
                tuple: (peak_width, peak_width_pixels, flux_range) values
            Raises:
                None
        """
        peak_width = PEAK_WIDTH_LOW_RES if self.resolution == 'LOW' else PEAK_WIDTH_HIGH_RES
        flux_range = 2 * peak_width
        angstroms_to_pixels = self.wavelength_data[1] - self.wavelength_data[0]
        peak_width_pixels = math.floor(peak_width / angstroms_to_pixels)

        return peak_width, peak_width_pixels, flux_range
    
    def smooth_data(self, sigma):
        """
            Applies Gaussian smoothing to the spectrum data.
            Arguments:
                sigma (float): Standard deviation for Gaussian kernel
            Returns:
                tuple: (smoothed_wavelength, smoothed_flux, smoothed_error) arrays
            Raises:
                None
        """
        smoothed_wavelength = gaussian_filter(self.wavelength_data, sigma)
        smoothed_flux = gaussian_filter(self.flux_data, sigma)
        smoothed_error = gaussian_filter(self.error_data, sigma)
        
        return smoothed_wavelength, smoothed_flux, smoothed_error
    
    def save_data(self, emission_lines):
        """
            Saves processed emission line data to various output formats.
            Arguments:
                emission_lines (EmissionLines): Object containing emission line data
            Returns:
                None
            Raises:
                None
        """
        data_list = self._prepare_data_list(emission_lines)
        data_header = self._prepare_data_header()
        
        self._save_ecsv(data_list, data_header)
        self._save_csv(data_list, data_header)
        self._save_fits(data_list, data_header)