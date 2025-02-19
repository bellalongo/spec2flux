from astropy.modeling.models import Voigt1D, Gaussian1D
from astropy.modeling import fitting, CompoundModel
from astropy.modeling.fitting import NonFiniteValueError
import numpy as np

from .model_fitter import ModelFitter, ModelFitterError

# ------------------------------
# Constants
# ------------------------------
TOLERANCE = 6  # Tolerance for grouping emission lines (in Angstroms)

# ------------------------------
# EmissionLines Class
# ------------------------------
class EmissionLines(object):
    def __init__(self, spectrum):
        """
        
        """
        self.spectrum = spectrum
        self.model_fitter = ModelFitter(spectrum)

        self.grouped_lines = self._group_emission_lines()

        # Initalize line list
        self.line_list = self._initialize_line_list()

        # Find the continuum
        self.spectrum.spectrum_continuum = self._find_spectrum_continuum()

        # Get emission line models
        self.line_models = self._get_line_models()

    # ------------------------------
    # Private Helper Functions
    # ------------------------------
    def _find_spectrum_continuum(self):
        """

        """
        # Check what kind !
        cont_arr = []

        # Use line_list instead of grouped_lines
        for ion in self.line_list:
            for emission_line in self.line_list[ion].values():
                group_mask = (
                    (self.spectrum.wavelength_data > emission_line.group_lam[0] - self.spectrum.line_width) &
                    (self.spectrum.wavelength_data < emission_line.group_lam[-1] + self.spectrum.line_width)
                )
                
                group_continuum = self._create_trendline(
                    self.spectrum.wavelength_data[group_mask],
                    self.spectrum.flux_data[group_mask]
                )
                cont_arr.extend(group_continuum)

        return np.percentile(cont_arr, 25)
    
    def _group_emission_lines(self):
        """

        """
        ion_groups = {}

        # Loop through emission lines
        for _, row in self.spectrum.rest_lam_data.iterrows():
            ion = row["Ion"]
            wavelength = float(row["Wavelength"])

            # Skip wavelengths outside the analysis range
            if wavelength < self.spectrum.min_wavelength:
                continue

            # Initialize ion group if it doesn't exist
            if ion not in ion_groups:
                ion_groups[ion] = [[wavelength]]
                continue

            # Check if the wavelength is close to any existing group
            close_group_found = False
            for group in ion_groups[ion]:
                if abs(max(group) - wavelength) <= TOLERANCE:
                    group.append(wavelength)
                    close_group_found = True
                    break

            # If no close group was found, create a new group
            if not close_group_found:
                ion_groups[ion].append([wavelength])

        return ion_groups

    def _create_trendline(self, wavelength_data, flux_data):
        """

        """
        length = len(wavelength_data) - 1
        flux_list_left = []
        flux_list_right = []

        # Collect flux values from the left and right of the peak
        for i in range(int(self.spectrum.line_width_pixels / 2)):
            flux_list_left.append(flux_data[i])
            flux_list_right.append(flux_data[length - i])

        # Calculate the average flux for the left and right
        avg_flux_left = np.mean(flux_list_left)
        avg_flux_right = np.mean(flux_list_right)

        # Use the lesser of the two values as the average flux
        avg_flux = min(avg_flux_left, avg_flux_right)

        # Create a continuum array
        continuum_array = [avg_flux for _ in range(length + 1)]
        return continuum_array
    
    def _initialize_line_list(self):
        """

        """
        line_list = {}

        # Iteratre through each ion
        for ion in self.grouped_lines:
            # Initialize an instance in linelist for the ion
            line_list[ion] = {}
            
            # Iterate through each group
            for group in self.grouped_lines[ion]:
                # Mask out of emission line group 
                group_mask = (self.spectrum.wavelength_data > group[0] - self.spectrum.line_width) & (
                    self.spectrum.wavelength_data < group[len(group) - 1] + self.spectrum.line_width) 
                
                # Check if data is valid
                if not any(self.spectrum.flux_data[group_mask]):
                    continue

                # Add emission line to line list
                emission_line = self.EmissionLine(
                    ion = ion, 
                    group_lam = group, 
                    obs_lam = None,
                    blended_bool = len(group) > 1,
                    doppler_candidate = False, 
                    noise_bool = True, # change name to is_noise
                    model_params = None,
                    continuum = None, 
                    flux_error = None
                )
                
                # Add to linelist
                line_list[ion][group[0]] = emission_line

        return line_list 
    
    def _get_line_models(self):
        """

        """
        line_models = {}

        for ion in self.line_list:
            for emission_line in self.line_list[ion].values():
                # Mask out the emission line group
                group_mask = (
                    (self.spectrum.wavelength_data > emission_line.group_lam[0] - self.spectrum.line_width) &
                    (self.spectrum.wavelength_data < emission_line.group_lam[-1] + self.spectrum.line_width)
                )

                if not any(self.spectrum.flux_data[group_mask]):
                    continue

                try:
                    # Use ModelFitting class instead of direct creation
                    compound_model = self.model_fitter.create_model(emission_line.group_lam)
                    
                    # Use ModelFitting class for fitting
                    fitted_model, fitter = self.model_fitter.fit_model(
                        compound_model,
                        self.spectrum.wavelength_data[group_mask],
                        self.spectrum.flux_data[group_mask]
                    )

                    # Use ModelFitting class to save parameters
                    model_params = self.model_fitter.save_model_params(fitted_model, fitter)
                    emission_line.model_params = model_params

                    # Store the fitted model and mask
                    line_models[tuple(emission_line.group_lam)] = {
                        'model': fitted_model,
                        'mask': group_mask,
                        'ion': ion
                    }

                except ModelFitterError:
                    continue

        return line_models
    
    # ------------------------------
    # Public Methods
    # ------------------------------
    def update_selected_lines(self, plotter, purpose):
        """
        
        """
        selected = {wavelength: data for wavelength, data in 
                   plotter.selected_lines.items() if data['selected']}
        
        for wavelength, line_data in selected.items():
            try:
                line = self.line_list[line_data['ion']][wavelength]
                if purpose == 'Doppler':
                    line.doppler_candidate = True
                    line.noise_bool = False
                elif purpose == 'Noise':
                    line.noise_bool = False
            except KeyError:
                continue

    def emission_line_to_dict(self, line_obj):
        """

        """
        return {
            "ion": line_obj.ion,
            "group_lam": line_obj.group_lam,
            "obs_lam": float(line_obj.obs_lam),
            "noise_bool": line_obj.noise_bool,
            "blended_bool": line_obj.blended_bool,
            "doppler_candidate": line_obj.doppler_candidate,
            "model_params": line_obj.model_params,
            "continuum": float(line_obj.continuum),
            "flux_error": line_obj.flux_error,
        }
    
    # ------------------------------
    # EmissionLine Inner Class
    # ------------------------------
    class EmissionLine:
        def __init__(self, ion, group_lam, obs_lam, noise_bool, blended_bool, doppler_candidate, model_params, continuum, flux_error):
            """
            Initialize the EmissionLine class.

            Parameters:
                ion (str): The ion name.
                group_lam (list): List of wavelengths in the group.
                obs_lam (float): Observed wavelength.
                noise_bool (bool): Whether the line is noise.
                blended_bool (bool): Whether the line is blended.
                doppler_candidate: Doppler shift candidate.
                model_params: Model parameters for fitting.
                continuum (float): Continuum value.
                flux_error (tuple): Flux and error values.
            """
            self.ion = ion
            self.group_lam = group_lam
            self.obs_lam = obs_lam
            self.noise_bool = noise_bool
            self.blended_bool = blended_bool
            self.doppler_candidate = doppler_candidate
            self.model_params = model_params
            self.continuum = continuum
            self.flux_error = flux_error