import numpy as np

# ------------------------------
# Constants
# ------------------------------
TOLERANCE = 10  # Tolerance for grouping emission lines (in Angstroms)

# ------------------------------
# EmissionLines Class
# ------------------------------
class EmissionLines(object):
    def __init__(self, spectrum):
        """
        
        """
        self.spectrum = spectrum
        self.grouped_lines = self._group_emission_lines()

        # Initialize emission line list
        self.line_list = []

        # Find the continuum
        self.spectrum_continuum = self._find_spectrum_continuum()

    # ------------------------------
    # Private Helper Functions
    # ------------------------------
    def _find_spectrum_continuum(self):
        """

        """
        cont_arr = []

        # Iterate through all of the rest wavelengths
        for ion in self.grouped_lines:
            for group in self.grouped_lines[ion]:
                # Mask out the grouped emission lines
                group_mask = (
                    (self.spectrum.wavelength_data > group[0] - self.spectrum.line_width) &
                    (self.spectrum.wavelength_data < group[-1] + self.spectrum.line_width)
                )

                # Find the continuum for each emission line
                group_continuum = self._create_trendline(
                    self.spectrum.wavelength_data[group_mask],
                    self.spectrum.flux_data[group_mask]
                )
                cont_arr += list(group_continuum)

        # Calculate the 25th percentile as the continuum
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
    
    # ------------------------------
    # Public Methods
    # ------------------------------
    def emission_line_to_dict(self, line_obj):
        """

        """
        return {
            "ion": line_obj.ion,
            "group_lam": line_obj.group_lam,
            "obs_lam": float(line_obj.obs_lam),
            "noise_bool": line_obj.noise_bool,
            "blended_bool": line_obj.blended_bool,
            "doppler_candidate": None,
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

    
    # def find_spectrum_continuum(self):
    #     """
            
    #     """
    #     cont_arr = []

    #     # Iterate through all of the rest wavelengths
    #     for ion in self.grouped_lines:
    #         for group in self.grouped_lines[ion]:
    #             # Mask out the grouped emission lines
    #             group_mask = (self.spectrum.wavelength_data > group[0] - self.spectrum.line_width) & (
    #                 self.spectrum.wavelength_data < group[len(group) - 1] + self.spectrum.line_width) 
                
    #             # Find the continuum for each emission line
    #             group_continuum = self.split_create_trendline(self.spectrum.wavelength_data[group_mask],
    #                                                           self.spectrum.flux_data[group_mask])
                
    #             cont_arr += list(group_continuum)

    #     # Calculate the mean to be the continuum
    #     return np.percentile(cont_arr, 25)



    # def grouping_emission_lines(self):
    #     """
    #         Groups blended emission lines together based off of ion, and a pre-determined tolerance
    #         (Note: this function assumes the DEM data is strictly increasing)
    #         Parameters: 
    #                     None
    #         Returns:
    #                     ion_groups: dictionary with ion name as the key, and  blended groups for that as the value
    #     """
    #     # Initialize variables
    #     tolerance = 1.5 # ADJUST ME (was 10 previously)
    #     ion_groups = {}
    #     close_group_found = False

    #     # Loop through emission lines
    #     for _, row in self.spectrum.rest_lam_data.iterrows():
    #         ion = row["Ion"]
    #         wavelength = float(row["Wavelength"])

    #         # Check if wavelength is within wavelength bounds
    #         if wavelength < self.spectrum.min_wavelength: 
    #             continue

    #         # Check if ion does not exist in the dictionary
    #         if ion not in ion_groups:
    #             ion_groups[ion] = [[wavelength]]
                
    #         # Ion is already in the dictionary
    #         else:
    #             close_group_found = False
    #             # Iterate through each group for that ion
    #             for group in ion_groups[ion]:
    #                 # If the largest value in the group - wavelength is less than the tolerance, add to that group
    #                 if abs(max(group) - wavelength) <= tolerance:
    #                     group.append(wavelength)
    #                     close_group_found = True
    #                     break
                
    #             # If no close group was found, add to that ion as a seperate group
    #             if not close_group_found:
    #                 ion_groups[ion].append([wavelength])

    #     return ion_groups
    

    # def split_create_trendline(self, wavelength_data, flux_data):
    #     """
    #         Create a continuum trendline using the average flux from the left and right of the peak
    #         Parameters: 
    #                     wavelength_data: masked wavelength data from the spectra     
    #                     flux_data: masked flux data from the spectra
    #         Returns:
    #                     continuum_array: continuum data for the current peak
    #     """
    #     # Initialize variables
    #     length = len(wavelength_data) - 1
    #     flux_list_left = []
    #     flux_list_right = []
        
    #     # Make an array of all flux that aren't included in the peak
    #     for i in range(0, int(self.spectrum.line_width_pixels/2)):
    #         flux_list_left.append(flux_data[i])
    #         flux_list_right.append(flux_data[length-i])
            
    #     # Find the average flux for the left and right
    #     avg_flux_left = sum(flux_list_left)/len(flux_list_left)
    #     avg_flux_right = sum(flux_list_right)/len(flux_list_right)
        
    #     # Use the lesser of the two values as the average flux
    #     if avg_flux_left < avg_flux_right:
    #         avg_flux = avg_flux_left
    #     else:
    #         avg_flux = avg_flux_right
        
    #     continuum_array = [avg_flux for _ in range(length + 1)]
        
    #     return continuum_array
    

    # def emission_line_to_dict(self, line_obj):
    #     """
    #         Saves an emission line object as a dictionary 
    #         Parameters: 
    #                     line_obj: emission line object
    #         Returns:
    #                     a dictionary of the emission line
    #     """
    #     return {
    #         "ion": line_obj.ion,
    #         "group_lam": line_obj.group_lam, 
    #         "obs_lam": float(line_obj.obs_lam),
    #         "noise_bool": line_obj.noise_bool,
    #         "blended_bool": line_obj.blended_bool,
    #         "doppler_candidate": None,
    #         "model_params": line_obj.model_params, # CHECK ME !!!!
    #         "continuum": float(line_obj.continuum),
    #         "flux_error": line_obj.flux_error
    #     }


    # class Emission_Line:
    #     def __init__(self, ion, group_lam, obs_lam, noise_bool, blended_bool, doppler_candidate, model_params, continuum, flux_error):
    #         self.ion = ion
    #         self.group_lam = group_lam
    #         self.obs_lam = obs_lam
    #         self.noise_bool = noise_bool
    #         self.blended_bool = blended_bool
    #         self.doppler_candidate = doppler_candidate
    #         self.model_params = model_params
    #         self.continuum = continuum
    #         self.flux_error = flux_error 