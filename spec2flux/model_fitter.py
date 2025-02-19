from astropy.modeling.models import Voigt1D, Gaussian1D
from astropy.modeling import fitting, CompoundModel
from astropy.modeling.fitting import NonFiniteValueError
import numpy as np

class ModelFitter:
    """

    """
    def __init__(self, spectrum):
        """

        """
        self.spectrum = spectrum

    # ------------------------------
    # Private Helper Methods
    # ------------------------------
    def _create_voigt_profile(self, line, line_mask):
        """
        
        """
        init_x0 = line
        init_amp = np.max(self.spectrum.flux_data[line_mask])
        init_fwhm_g = self.spectrum.line_width/5
        init_fwhm_l = self.spectrum.line_width/5

        return Voigt1D(
            x_0=init_x0,
            amplitude_L=init_amp,
            fwhm_L=init_fwhm_l,
            fwhm_G=init_fwhm_g
        )

    def _create_gaussian_profile(self, line, line_mask):
        """

        """
        init_amp = np.max(self.spectrum.flux_data[line_mask])
        init_mean = line
        init_stddev = self.spectrum.line_width / 10

        return Gaussian1D(
            amplitude=init_amp,
            mean=init_mean,
            stddev=init_stddev
        )

    def _combine_models(self, model_profiles):
        """

        """
        if not model_profiles:
            raise ValueError("No model profiles to combine")
            
        compound_model = model_profiles[0]
        for profile in model_profiles[1:]:
            compound_model += profile
        return compound_model

    def _extract_voigt_params(self, fitted_model, param_errors, has_errors):
        """

        """
        params = []
        
        if isinstance(fitted_model, CompoundModel):
            components = fitted_model
        else:
            components = [fitted_model]

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

        """
        params = []
        
        if isinstance(fitted_model, CompoundModel):
            components = fitted_model
        else:
            components = [fitted_model]

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

        """
        model_profiles = []

        for line in group:
            # Mask out each emission line
            line_mask = (self.spectrum.wavelength_data > line - self.spectrum.line_width/2) & (
                self.spectrum.wavelength_data < line + self.spectrum.line_width/2)

            if self.spectrum.line_fit_model == 'Voigt':
                profile = self._create_voigt_profile(line, line_mask)
            else:
                profile = self._create_gaussian_profile(line, line_mask)
                
            model_profiles.append(profile)

        return self._combine_models(model_profiles)
    
    def create_model_profile(self, model_params):
        """

        """
        model_profiles = []

        for params in model_params:
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

        """
        try:
            fitter = fitting.LevMarLSQFitter()
            fitted_model = fitter(model, wavelength_data, flux_data)
            return fitted_model, fitter
        except (RuntimeError, TypeError, NonFiniteValueError) as e:
            raise ModelFitterError(f"Failed to fit model: {str(e)}")

    def save_model_params(self, fitted_model, fitter):
        """

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

        if self.spectrum.line_fit_model == 'Voigt':
            # Check model type
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

        else:  # Gaussian
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

    """
    pass