import astropy.units as u
from astropy.modeling.models import Voigt1D
from astropy.modeling import fitting, CompoundModel
from astropy.modeling.fitting import NonFiniteValueError
import json
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import exists
import seaborn as sns
import sys


class FluxCalculator(object): 
    def __init__(self, spectrum, emission_lines, fresh_start):
        self.spectrum = spectrum
        self.emission_lines = emission_lines
        self.fresh_start = fresh_start

        # Initialize boolean lists to be used in flux calculation
        self.doppler_bool_list = []
        self.noise_bool_list = []

        self.flux_calc()

    
    def flux_calc(self):
        '''
        '''
        # Check run type
        if not self.fresh_start:
            # Load saved doppler shift
            doppler_shift = np.loadtxt(self.spectrum.doppler_dir)*(u.km/u.s)

            # Read emission_line data from JSON file
            with open(self.spectrum.emission_lines_dir, "r") as json_file:
                emission_line_data = json.load(json_file)

            # Reconstruct emission_line objects from dictionaries
            self.emission_lines.line_list.clear() 
            for data in emission_line_data:
                emission_line_obj = self.emission_lines.Emission_Line(**data) # Check if this is legal -> it isn't fix me 
                self.emission_lines.line_list.append(emission_line_obj)
        else:
            if exists(self.spectrum.doppler_dir):
                # Remove all files
                os.remove(self.spectrum.doppler_dir)
                os.remove(self.spectrum)
                os.remove(self.spectrum.fits_dir)
                os.remove(self.spectrum.ecsv_dir)
                os.remove(self.spectrum.csv_dir)
                os.remove(self.spectrum.final_plot_dir)
            # else: # MAYBE ADD ME ??? WHEN DONE !
            #     print('Either try to run again, or change run type.')

            # Calculate doppler
            doppler_shift = self.doppler_shift_calc()

    
    def doppler_selection_plots(self):
        # Store line_list index of doppler candidates
        self.candidate_indices = [] # MAYBE delete

        # Iterate through groups
        for ion in self.emission_lines.grouped_lines:
            for group in self.emission_lines.grouped_lines[ion]:
                # Store voigt profiles for each group
                voigt_profiles = []

                # Mask out the emission line group
                group_mask = (self.spectrum.wavelength_data > group[0] - self.spectrum.line_width) & (
                    self.spectrum.wavelength_data < group[len(group) - 1] + self.spectrum.line_width) 

                # Check if data is valid
                if not any(self.spectrum.flux_data[group_mask]):
                    continue

                # Iterate through each emission line in the group
                for line in group:
                    # Mask out each emission line
                    line_mask = (self.spectrum.wavelength_data > line - self.spectrum.line_width/2) & (
                        self.spectrum.wavelength_data < line + self.spectrum.line_width/2) 
                    
                    # Voigt initial parama guesses
                    init_x0 = line
                    init_amp = np.max(self.spectrum.flux_data[line_mask]) 
                    init_fwhm_g = self.spectrum.line_width/5
                    init_fwhm_l = self.spectrum.line_width/5

                    # Voigt profile
                    voigt_profile = Voigt1D(x_0 = init_x0, amplitude_L = init_amp, fwhm_L = init_fwhm_l, fwhm_G = init_fwhm_g)
                    voigt_profiles.append(voigt_profile)
                
                # Combine voigt profiles
                composite_model = voigt_profiles[0]
                for voigt_profile in voigt_profiles[1:]:
                    composite_model += voigt_profile

                # Update emission line list
                line_group = self.emission_lines.Emission_Line(
                    ion = ion, 
                    group_lam = group, 
                    obs_lam = None, 
                    noise_bool = None, 
                    blended_bool = len(group) > 1, 
                    doppler_candidate = None, 
                    model_params = None,
                    continuum = None,
                    flux_error = None
                    )
                self.emission_lines.line_list.append(line_group)

                # Try to fit the model
                try: 
                    fitter = fitting.LevMarLSQFitter()
                    fitted_model = fitter(composite_model, self.spectrum.wavelength_data[group_mask], 
                                          self.spectrum.flux_data[group_mask])
                except (RuntimeError, TypeError, NonFiniteValueError):
                    continue

                # Save model info
                curr_index = len(self.emission_lines.line_list) - 1
                self.candidate_indices.append(curr_index) # maybe just check if there is a fitted_model instead ??

                # Save fitted model params
                model_params = []
                
                if isinstance(fitted_model, CompoundModel):
                    for component in fitted_model:
                        model_params.append({
                            'x_0': component.x_0.value,
                            'amplitude_L': component.amplitude_L.value,
                            'fwhm_L': component.fwhm_L.value,
                            'fwhm_G': component.fwhm_G.value
                        })
                else:
                    model_params.append({
                        'x_0': fitted_model.x_0.value,
                        'amplitude_L': fitted_model.amplitude_L.value,
                        'fwhm_L': fitted_model.fwhm_L.value,
                        'fwhm_G': fitted_model.fwhm_G.value
                    })
                                            
                # Save model params
                self.emission_lines.line_list[curr_index].model_params = model_params # MAYBE save params instead?

                # Plot basics
                sns.set_theme()
                fig = plt.figure(figsize=(14,7), facecolor="white")
                ax = fig.add_subplot()
                plt.title(f"Flux vs Wavelength for {ion}")
                fig.suptitle("Click 'y' if should be used for doppler calculation, 'n' if not", fontweight='bold')
                plt.xlabel('Wavelength (Å)', fontsize =12)
                plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=12)
                cid = fig.canvas.mpl_connect('key_press_event', lambda event: self.on_key(event, 'Doppler Calculation'))

                plt.plot(self.spectrum.wavelength_data[group_mask], self.spectrum.flux_data[group_mask], linewidth=1)
                voigt_fit, = plt.plot(self.spectrum.wavelength_data[group_mask], fitted_model(self.spectrum.wavelength_data[group_mask]), color = "#111D4A")     

                # Plotting rest and obs wavelengths
                for i, wavelength in enumerate(group):
                    # Rest and observed wavlengths
                    rest_lam = plt.axvline(x = wavelength, color = "#F96E46", linestyle=((0, (5, 5))), linewidth=1)
                    obs_lam = plt.axvline(x=fitted_model.x_0.value if len(group) == 1 else fitted_model[i].x_0.value,
                                        color="#8F1423", linewidth=1)
                    
                legend = plt.legend([rest_lam, obs_lam, voigt_fit], ["Rest Wavelength", "Observed Wavelength", "Voigt Fit"])
                legend.get_frame().set_facecolor('white')
                plt.show()

        assert len(self.doppler_bool_list) > 0, "You didn't click 'y' and 'n' to choose did you?"
        assert len(self.doppler_bool_list) == len(self.candidate_indices), "Did you click out of the figure instead of clicking 'y' or 'n'?" 

        # # Calculate doppler shift
        # dv = []
        # # Iterate through each doppler candidate
        # for i, boolean in enumerate(doppler_bool_list):
        #     group_doppler = []
            
        #     # If is a doppler candidate, set as one in dict
        #     if boolean: 
        #         emission_line_objs[i].update_doppler_candidate(True)

        #         # Calculate doppler for each emission line in that group
        #         for j, rest_wavelength in enumerate(rest_candidates[i]):
        #             u_rest_lam = rest_wavelength * u.AA
        #             u_obs_lam = obs_candidates[i][j] * u.AA
        #             group_doppler.append(u_obs_lam.to(u.km/u.s,  equivalencies=u.doppler_optical(u_rest_lam)))
        #         dv.append(sum(group_doppler)/ len(group_doppler))

        # doppler_shift = sum(dv)/len(dv)

        # # Store value
        # with open(doppler_filename, 'a') as f:
        #     f.write(f"{doppler_shift.value:.3f}\n")
        
        # return doppler_shift

    def doppler_shift_calc(self):
        # Show doppler shift plots
        self.doppler_selection_plots()

        # 
    

    def on_key(self, event, purpose):
        """
            Event function that determines if a key was clicked
            Parameters: 
                        event: key press event
                        purpose: either doppler or noise
            Returns:
                        None
        """
        valid_keys = {'y', 'n'}

        if purpose not in ["Noise Detection", "Doppler Calculation"]:
            sys.exit("Invalid purpose, select 'Noise Detection' or 'Doppler Calculation'")

        if event.key not in valid_keys:
            print("Invalid key input, select 'y' or 'n'")

        if purpose == "Noise Detection":
            self.noise_bool_list.append(event.key == 'y')
            plt.close()

        elif purpose == "Doppler Calculation":
            self.doppler_bool_list.append(event.key == 'y')
            plt.close()



    