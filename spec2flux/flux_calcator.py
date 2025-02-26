import astropy.units as u
import json
import numpy as np
import os

from .model_fitter import ModelFitter

# ------------------------------
# Constants
# ------------------------------
UPPER_LIMIT_FACTOR = 3  # Factor for calculating upper limit from error

# ------------------------------
# FluxCalculator Class
# ------------------------------
class FluxCalculator:
    """
        This class calculates fluxes for emission lines in a spectrum.
        It handles both modeled lines and unmodeled lines, applies continuum subtraction,
        and provides methods to save results in various formats.
        Attributes:
            spectrum (SpectrumData): The spectrum data object
            emission_lines (EmissionLines): Object containing emission line data
            plotter (Plotter): Object for plotting visualization
            model_fitter (ModelFitter): Object for model fitting
    """
    
    def __init__(self, spectrum, emission_lines):
        """
            Initializes the FluxCalculator with spectrum and emission line data.
            Arguments:
                spectrum (SpectrumData): The spectrum data object
                emission_lines (EmissionLines): Object containing emission line data
                plotter (Plotter): Object for plotting visualization
        """
        self.spectrum = spectrum
        self.emission_lines = emission_lines
        self.model_fitter = ModelFitter(spectrum)

    # ------------------------------
    # Private Helper Functions
    # ------------------------------
    def _calculate_continuum(self, line, wavelength_data, flux_data):
        """
            Calculates the continuum level for a specific emission line.
            Arguments:
                line (EmissionLine): The emission line object
                wavelength_data (ndarray): Array of wavelength values
                flux_data (ndarray): Array of flux values
            Returns:
                list: Array of continuum values matching the input wavelength array length
        """
        # Check if fitting the entire spectrum as the continuum
        if self.spectrum.cont_fit == 'Complete':
            continuum = [self.spectrum.spectrum_continuum] * len(wavelength_data)

        # Individual continuum
        else:
            # Check if there was a successful fit
            if line.model_params:
                # Use fit min as continuum
                model_profile = self.model_fitter.create_model_profile(line.model_params)
                continuum = [min(model_profile(wavelength_data))] * len(wavelength_data)
            else:
                # Use trendline as continuum
                continuum = self.emission_lines._create_trendline(wavelength_data, flux_data)

        line.continuum = continuum[0]

        return continuum

    def _calculate_wavelength_edges(self, group_wavelength_data):
        """
            Calculates the wavelength bin edges for numerical integration.
            Arguments:
                group_wavelength_data (ndarray): Array of wavelength values
            Returns:
                tuple: (w0, w1) Arrays of lower and upper wavelength bin edges
        """
        # Calculate difference between each wavelength point
        diff = np.diff(group_wavelength_data)

        # Make a boxes correlating to these edges 
        diff0 = np.concatenate((np.array([diff[0]]), diff))
        diff1 = np.concatenate((diff, np.array([diff[-1]])))

        w0 = group_wavelength_data - diff0/2.
        w1 = group_wavelength_data + diff1/2.
        
        return w0, w1

    def _calculate_single_flux(self, emission_line):
        """
            Calculates the flux for a single emission line or line group.
            Arguments:
                emission_line (EmissionLine): The emission line object
            Returns:
                tuple: (flux, error) The calculated flux and its associated error
        """
        # Create mask for the line
        group_mask = (
            (self.spectrum.wavelength_data > emission_line.group_lam[0] - self.spectrum.line_width) &
            (self.spectrum.wavelength_data < emission_line.group_lam[-1] + self.spectrum.line_width)
        )
        
        # Get masked data
        wavelength_data = self.spectrum.wavelength_data[group_mask]
        flux_data = self.spectrum.flux_data[group_mask]
        error_data = self.spectrum.error_data[group_mask]

        # Update line information based on existing doppler shift
        rest_lam = emission_line.group_lam[-1] * u.AA
        emission_line.obs_lam = self.spectrum.doppler_shift.to(
            u.AA, equivalencies=u.doppler_optical(rest_lam)).value

        # Calculate continuum and edges
        w0, w1 = self._calculate_wavelength_edges(wavelength_data)
        continuum = self._calculate_continuum(emission_line, wavelength_data, flux_data)
        
        # Calculate flux components
        if emission_line.model_params:
            model_profile = self.model_fitter.create_model_profile(emission_line.model_params)
            total_sumflux = np.sum((model_profile(wavelength_data))*(w1 - w0))
        else:
            total_sumflux = np.sum(flux_data*(w1-w0))

        # Calculate final flux and error
        continuum_sumflux = np.sum(continuum*(w1 - w0))
        sumerror = (np.sum(error_data**2 * (w1 - w0)**2))**0.5
        
        if emission_line.noise_bool:
            return -UPPER_LIMIT_FACTOR * sumerror, 0
        else:
            return total_sumflux - continuum_sumflux, sumerror

    def _save_flux_results(self, line_dicts):
        """
            Saves emission line data to a JSON file.
            Arguments:
                line_dicts (list): List of emission line dictionaries
            Returns:
                None
        """
        with open(self.spectrum.emission_lines_dir, "w") as json_file:
            json.dump(line_dicts, json_file, indent=4)

    def _delete_existing_files(self):
        """
            Deletes existing output files to ensure clean results.
            Arguments:
                None
            Returns:
                None
        """
        files_to_delete = [
            self.spectrum.emission_lines_dir,
            self.spectrum.fits_dir,
            self.spectrum.ecsv_dir,
            self.spectrum.csv_dir,
            self.spectrum.final_plot_dir
        ]
        
        for file_path in files_to_delete:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
            except OSError as e:
                print(f"Warning: Could not delete {file_path}: {e}")

    # ------------------------------
    # Public Methods
    # ------------------------------
    def process_spectrum(self):
        """
            Processes the spectrum by calculating fluxes for all emission lines and saving results.
            Arguments:
                None
            Returns:
                None
        """
        line_dicts = []

        # Iterate through each emission line group
        for ion in self.emission_lines.line_list:
            for emission_line in self.emission_lines.line_list[ion].values():
                # Grab flux and error
                flux, error = self._calculate_single_flux(emission_line)
                emission_line.flux_error = (flux, error)

                # Append emission line to line_dicts
                line_dicts.append(self.emission_lines.emission_line_to_dict(emission_line))

        # Save emission line information
        self._save_flux_results(line_dicts)
        self.spectrum.save_data(self.emission_lines)