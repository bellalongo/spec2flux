import astropy.units as u
from astropy.constants import c
import numpy as np

class DopplerCalculator:
    """
        This class calculates the Doppler shift from emission lines in a spectrum.
        It uses the model parameters from fitted emission lines to determine individual line velocities,
        then computes a weighted mean to obtain the final Doppler shift value.
        Attributes:
            spectrum (SpectrumData): The spectrum data object
            doppler_shift (Quantity): The calculated Doppler shift in km/s
            doppler_error (Quantity): The error in the calculated Doppler shift in km/s
    """
    def __init__(self, spectrum):
        """
            Initializes the DopplerCalculator with a spectrum object.
            Arguments:
                spectrum (SpectrumData): The spectrum data object containing wavelength and flux data
        """
        self.spectrum = spectrum
        self.doppler_shift = None
        self.doppler_error = None

    # ------------------------------
    # Private Helper Methods
    # ------------------------------
    def _compute_weighted_mean(self, values, errors):
        """
            Computes the weighted mean of a set of values with their associated errors.
            Arguments:
                values (list): List of values as Astropy Quantity objects
                errors (list): List of error values as Astropy Quantity objects
            Returns:
                tuple: (weighted_mean, weighted_error) as Astropy Quantity objects in km/s
        """
        if len(values) > 1:
            weights = [1/error.value**2 for error in errors]
            weighted_mean = sum(v.value * w for v, w in zip(values, weights)) / sum(weights)
            weighted_error = np.sqrt(1 / sum(weights))
        else:
            weighted_mean = values[0].value
            weighted_error = errors[0].value

        return weighted_mean * u.km/u.s, weighted_error * u.km/u.s

    def _save_results(self):
        """
            Saves the calculated Doppler shift and error to a file.
            Arguments:
                None
            Returns:
                None
        """
        with open(self.spectrum.doppler_dir, 'a') as f:
            f.write(
                f"{self.doppler_shift.value:.3f} ± {self.doppler_error.value:.3f}\n"
            )

    def _calculate_single_line_velocity(self, rest_wavelength, model_params):
        """
            Calculates the velocity for a single emission line.
            Arguments:
                rest_wavelength (float): Rest wavelength of the emission line in Angstroms
                model_params (dict): Model parameters from the fitted line profile
            Returns:
                tuple: (velocity, velocity_error) as Astropy Quantity objects in km/s
        """
        rest_lam = rest_wavelength * u.AA
        
        # Get observed wavelength from model parameters
        if self.spectrum.line_fit_model == 'Voigt':
            obs_lam = model_params['x_0'] * u.AA
            obs_error = model_params['x_0_error'] * u.AA
        else:
            obs_lam = model_params['mean'] * u.AA
            obs_error = model_params['mean_error'] * u.AA

        # Calculate velocity
        velocity = obs_lam.to(
            u.km/u.s, 
            equivalencies=u.doppler_optical(rest_lam)
        )
        
        # Calculate velocity error
        velocity_error = abs(
            c.to(u.km/u.s).value * obs_error.value / rest_wavelength
        ) * u.km/u.s

        return velocity, velocity_error

    # ------------------------------
    # Public Methods
    # ------------------------------
    def calculate_doppler_shift(self, emission_lines):
        """
            Calculates the Doppler shift based on selected emission lines.
            Arguments:
                emission_lines (EmissionLines): Object containing emission line data
            Returns:
                Quantity: The calculated Doppler shift as an Astropy Quantity in km/s
            Raises:
                DopplerCalculationError: If no lines are selected for Doppler calculation
        """
        velocities = []
        velocity_errors = []

        for ion in emission_lines.line_list:
            for line in emission_lines.line_list[ion].values():
                if line.doppler_candidate:
                    line_velocity, line_error = self.calculate_line_velocity(line)
                    velocities.append(line_velocity)
                    velocity_errors.append(line_error)

        if not velocities:
            raise DopplerCalculationError("No lines selected for Doppler calculation")

        self.doppler_shift, self.doppler_error = self._compute_weighted_mean(
            velocities, velocity_errors
        )
        
        self._save_results()
        return self.doppler_shift

    def calculate_line_velocity(self, emission_line):
        """
            Calculates the velocity for a group of emission lines (possibly blended).
            Arguments:
                emission_line (EmissionLine): Emission line object
            Returns:
                tuple: (velocity, velocity_error) as Astropy Quantity objects in km/s
        """
        line_velocities = []
        line_errors = []

        for i, rest_wavelength in enumerate(emission_line.group_lam):
            velocity, velocity_error = self._calculate_single_line_velocity(
                rest_wavelength, 
                emission_line.model_params[i]
            )
            line_velocities.append(velocity)
            line_errors.append(velocity_error)

        # Compute weighted mean for the line group
        return self._compute_weighted_mean(line_velocities, line_errors)
    
class DopplerCalculationError(Exception):
    """
        Exception raised when there is an error in Doppler shift calculation.
        ** This typically occurs when no emission lines are selected for the calculation.
    """
    pass