from astropy.modeling.models import Voigt1D, Gaussian1D
from astropy.modeling import fitting, CompoundModel
from astropy.modeling.fitting import NonFiniteValueError
import numpy as np

class ModelFitter:
    """
        This class handles the creation and fitting of spectral line models to observed emission lines.
        It supports both Voigt and Gaussian profile models, and manages the fitting process and
        parameter extraction.
        Attributes:
            spectrum (SpectrumData): The spectrum data object containing configuration and data
    """
    def __init__(self, spectrum):
        """
            Initializes the ModelFitter with a spectrum object.
            Arguments:
                spectrum (SpectrumData): The spectrum data object
        """
        self.spectrum = spectrum

    # ------------------------------
    # Private Helper Methods
    # ------------------------------
    def _create_voigt_profile(self, line, line_mask):
        """
            Creates a Voigt profile model for a single emission line.
            Arguments:
                line (float): Wavelength of the emission line
                line_mask (ndarray): Boolean mask for the wavelength range of the line
            Returns:
                Voigt1D: Voigt profile model with initial parameters
        """
        # Create voigt profile initalial guesses 
        init_x0 = line
        init_amp = np.max(self.spectrum.flux_data[line_mask])
        init_fwhm_g = self.spectrum.line_width/5
        init_fwhm_l = self.spectrum.line_width/5

        # Create a voigt profile using initalized values 
        return Voigt1D(
            x_0=init_x0,
            amplitude_L=init_amp,
            fwhm_L=init_fwhm_l,
            fwhm_G=init_fwhm_g
        )

    def _create_gaussian_profile(self, line, line_mask):
        """
            Creates a Gaussian profile model for a single emission line.
            Arguments:
                line (float): Wavelength of the emission line
                line_mask (ndarray): Boolean mask for the wavelength range of the line
            Returns:
                Gaussian1D: Gaussian profile model with initial parameters
        """
        # Create gaussian profile initalial guesses 
        init_amp = np.max(self.spectrum.flux_data[line_mask])
        init_mean = line
        init_stddev = self.spectrum.line_width / 10

        # Create a Gaussian profile using initialized values 
        return Gaussian1D(
            amplitude=init_amp,
            mean=init_mean,
            stddev=init_stddev
        )

    def _combine_models(self, model_profiles):
        """
            Combines multiple model profiles into a compound model.
            Arguments:
                model_profiles (list): List of model profile objects
            Returns:
                CompoundModel: Combined model or single model if only one profile
            Raises:
                ValueError: If model_profiles is empty
        """
        if not model_profiles:
            raise ValueError("No model profiles to combine")
            
        # Combine the model profiles
        compound_model = model_profiles[0]
        for profile in model_profiles[1:]:
            compound_model += profile

        return compound_model

    def _extract_voigt_params(self, fitted_model, param_errors, has_errors):
        """
            Extracts parameters from a fitted Voigt model.
            Arguments:
                fitted_model: The fitted model (can be CompoundModel or single model)
                param_errors (ndarray): Array of parameter errors
                has_errors (bool): Whether error estimates are available
            Returns:
                list: List of dictionaries containing model parameters and their errors
        """
        params = []
        
        # Check if compound model
        if isinstance(fitted_model, CompoundModel):
            components = fitted_model
        else:
            components = [fitted_model]

        # Iterate through each component in the models and extract params
        for i, component in enumerate(components):
            params.append({
                'x_0': component.x_0.value,
                'x_0_error': param_errors[i*4] if has_errors else 0.1,
                'amplitude_L': component.amplitude_L.value,
                'amplitude_L_error': param_errors[i*4 + 1] if has_errors else 0.1,
                'fwhm_L': component.fwhm_L.value,
                'fwhm_L_error': param_errors[i*4 + 2] if has_errors else 0.1,
                'fwhm_G': component.fwhm_G.value,
                'fwhm_G_error': param_errors[i*4 + 3] if has_errors else 0.1
            })
        
        return params

    def _extract_gaussian_params(self, fitted_model, param_errors, has_errors):
        """
            Extracts parameters from a fitted Gaussian model.
            Arguments:
                fitted_model: The fitted model (can be CompoundModel or single model)
                param_errors (ndarray): Array of parameter errors
                has_errors (bool): Whether error estimates are available
            Returns:
                list: List of dictionaries containing model parameters and their errors
        """
        params = []
        
        # Check if compound model
        if isinstance(fitted_model, CompoundModel):
            components = fitted_model
        else:
            components = [fitted_model]

        # Iterate through each component in the models and extract params
        for i, component in enumerate(components):
            params.append({
                'amplitude': component.amplitude.value,
                'amplitude_error': param_errors[i*3] if has_errors else 0.1,
                'mean': component.mean.value,
                'mean_error': param_errors[i*3 + 1] if has_errors else 0.1,
                'stddev': component.stddev.value,
                'stddev_error': param_errors[i*3 + 2] if has_errors else 0.1
            })
        
        return params

    # ------------------------------
    # Public Methods
    # ------------------------------
    def create_model(self, group):
        """
            Creates a model for a group of emission lines.
            Arguments:
                group (list): List of wavelengths for the emission line group
            Returns:
                CompoundModel: Combined model for all lines in the group
            Raises:
                ValueError: If no model profiles could be created
        """
        model_profiles = []

        # Iterate through each line in the group
        for line in group:
            # Mask out each emission line
            line_mask = (self.spectrum.wavelength_data > line - self.spectrum.line_width/2) & (
                self.spectrum.wavelength_data < line + self.spectrum.line_width/2)

            # Check model type
            if self.spectrum.line_fit_model == 'Voigt':
                profile = self._create_voigt_profile(line, line_mask)
            else:
                profile = self._create_gaussian_profile(line, line_mask)
                
            model_profiles.append(profile)

        return self._combine_models(model_profiles)
    
    def create_model_profile(self, model_params):
        """
            Creates a model profile from saved parameters.
            Arguments:
                model_params (list): List of dictionaries containing model parameters
            Returns:
                CompoundModel: Combined model constructed from the parameters
            Raises:
                ValueError: If no model profiles could be created
        """
        model_profiles = []

        # Iterate through model parameters 
        for params in model_params:
            # Check model type
            if self.spectrum.line_fit_model == 'Voigt':
                model_profile = Voigt1D(
                    x_0=params['x_0'],
                    amplitude_L=params['amplitude_L'],
                    fwhm_L=params['fwhm_L'],
                    fwhm_G=params['fwhm_G']
                )
            else:
                model_profile = Gaussian1D(
                    amplitude=params['amplitude'],
                    mean=params['mean'],
                    stddev=params['stddev']
                )
            model_profiles.append(model_profile)

        return self._combine_models(model_profiles)

    def fit_model(self, model, wavelength_data, flux_data):
        """
            Fits a model to the observed data.
            Arguments:
                model: The model to fit (CompoundModel or single model)
                wavelength_data (ndarray): Array of wavelength values
                flux_data (ndarray): Array of flux values
            Returns:
                tuple: (fitted_model, fitter) The fitted model and the fitter object
            Raises:
                ModelFitterError: If the fitting process fails
        """
        try:
            fitter = fitting.LevMarLSQFitter()
            fitted_model = fitter(model, wavelength_data, flux_data)
            return fitted_model, fitter
        except (RuntimeError, TypeError, NonFiniteValueError) as e:
            raise ModelFitterError(f"Failed to fit model: {str(e)}")

    def save_model_params(self, fitted_model, fitter):
        """
            Extracts and saves parameters from a fitted model.
            Arguments:
                fitted_model: The fitted model (CompoundModel or single model)
                fitter: The fitter object used for fitting
            Returns:
                list: List of dictionaries containing model parameters and their errors
        """
        model_params = []
        default_error = 0.1

        # Check if we have valid parameter covariance information
        try:
            has_errors = (hasattr(fitter, 'fit_info') and 
                        'param_cov' in fitter.fit_info and 
                        fitter.fit_info['param_cov'] is not None)
            if has_errors:
                param_errors = np.sqrt(np.diag(fitter.fit_info['param_cov']))
            else:
                param_errors = None
        except (ValueError, TypeError, np.linalg.LinAlgError):
            has_errors = False
            param_errors = None

        # Voigt model
        if self.spectrum.line_fit_model == 'Voigt':
            # Check if compound model
            if isinstance(fitted_model, CompoundModel):
                for i, component in enumerate(fitted_model):
                    base_idx = i * 4  # 4 parameters per Voigt component
                    model_params.append({
                        'x_0': component.x_0.value,
                        'x_0_error': param_errors[base_idx] if param_errors is not None else default_error,
                        'amplitude_L': component.amplitude_L.value,
                        'amplitude_L_error': param_errors[base_idx + 1] if param_errors is not None else default_error,
                        'fwhm_L': component.fwhm_L.value,
                        'fwhm_L_error': param_errors[base_idx + 2] if param_errors is not None else default_error,
                        'fwhm_G': component.fwhm_G.value,
                        'fwhm_G_error': param_errors[base_idx + 3] if param_errors is not None else default_error
                    })
            else:
                model_params.append({
                    'x_0': fitted_model.x_0.value,
                    'x_0_error': param_errors[0] if param_errors is not None else default_error,
                    'amplitude_L': fitted_model.amplitude_L.value,
                    'amplitude_L_error': param_errors[1] if param_errors is not None else default_error,
                    'fwhm_L': fitted_model.fwhm_L.value,
                    'fwhm_L_error': param_errors[2] if param_errors is not None else default_error,
                    'fwhm_G': fitted_model.fwhm_G.value,
                    'fwhm_G_error': param_errors[3] if param_errors is not None else default_error
                })

        # Gaussian model 
        else:  
            # Check if compound model
            if isinstance(fitted_model, CompoundModel):
                for i, component in enumerate(fitted_model):
                    base_idx = i * 3  # 3 parameters per Gaussian component
                    model_params.append({
                        'amplitude': component.amplitude.value,
                        'amplitude_error': param_errors[base_idx] if param_errors is not None else default_error,
                        'mean': component.mean.value,
                        'mean_error': param_errors[base_idx + 1] if param_errors is not None else default_error,
                        'stddev': component.stddev.value,
                        'stddev_error': param_errors[base_idx + 2] if param_errors is not None else default_error
                    })
            else:
                model_params.append({
                    'amplitude': fitted_model.amplitude.value,
                    'amplitude_error': param_errors[0] if param_errors is not None else default_error,
                    'mean': fitted_model.mean.value,
                    'mean_error': param_errors[1] if param_errors is not None else default_error,
                    'stddev': fitted_model.stddev.value,
                    'stddev_error': param_errors[2] if param_errors is not None else default_error
                })

        return model_params


class ModelFitterError(Exception):
    """
        Exception raised when there is an error in the model fitting process.
        ** This can occur due to convergence issues, numerical errors, or invalid input data.
    """
    pass