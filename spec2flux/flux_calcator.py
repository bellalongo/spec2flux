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

    """
    
    def __init__(self, spectrum, emission_lines, plotter):
        """

        """
        self._spectrum = spectrum
        self._emission_lines = emission_lines
        self._plotter = plotter
        self._model_fitter = ModelFitter(spectrum)

    # ------------------------------
    # Private Helper Functions
    # ------------------------------
    def _calculate_continuum(self, line, wavelength_data, flux_data):
        """

        """
        if self._spectrum.cont_fit == 'Complete':
            continuum = [self._spectrum.spectrum_continuum] * len(wavelength_data)
        else:
            if line.model_params:
                model_profile = self._model_fitter.create_model_profile(line.model_params)
                continuum = [min(model_profile(wavelength_data))] * len(wavelength_data)
            else:
                continuum = self._emission_lines._create_trendline(wavelength_data, flux_data)

        line.continuum = continuum[0]
        return continuum

    def _calculate_wavelength_edges(self, group_wavelength_data):
        """

        """
        diff = np.diff(group_wavelength_data)
        diff0 = np.concatenate((np.array([diff[0]]), diff))
        diff1 = np.concatenate((diff, np.array([diff[-1]])))
        w0 = group_wavelength_data - diff0/2.
        w1 = group_wavelength_data + diff1/2.
        return w0, w1

    def _calculate_single_flux(self, emission_line):
        """

        """
        # Create mask for the line
        group_mask = (
            (self._spectrum.wavelength_data > emission_line.group_lam[0] - self._spectrum.line_width) &
            (self._spectrum.wavelength_data < emission_line.group_lam[-1] + self._spectrum.line_width)
        )
        
        # Get masked data
        wavelength_data = self._spectrum.wavelength_data[group_mask]
        flux_data = self._spectrum.flux_data[group_mask]
        error_data = self._spectrum.error_data[group_mask]

        # Update line information based on existing doppler shift
        rest_lam = emission_line.group_lam[-1] * u.AA
        emission_line.obs_lam = self._spectrum.doppler_shift.to(
            u.AA, equivalencies=u.doppler_optical(rest_lam)).value

        # Calculate continuum and edges
        w0, w1 = self._calculate_wavelength_edges(wavelength_data)
        continuum = self._calculate_continuum(emission_line, wavelength_data, flux_data)
        
        # Calculate flux components
        if emission_line.model_params:
            model_profile = self._model_fitter.create_model_profile(emission_line.model_params)
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

        """
        with open(self._spectrum.emission_lines_dir, "w") as json_file:
            json.dump(line_dicts, json_file, indent=4)

    def _delete_existing_files(self):
        """

        """
        files_to_delete = [
            self._spectrum.emission_lines_dir,
            self._spectrum.fits_dir,
            self._spectrum.ecsv_dir,
            self._spectrum.csv_dir,
            self._spectrum.final_plot_dir
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
            
        """
        line_dicts = []

        for ion in self._emission_lines.line_list:
            for emission_line in self._emission_lines.line_list[ion].values():
                flux, error = self._calculate_single_flux(emission_line)
                emission_line.flux_error = (flux, error)
                line_dicts.append(self._emission_lines.emission_line_to_dict(emission_line))

        self._save_flux_results(line_dicts)
        self._spectrum.save_data(self._emission_lines)