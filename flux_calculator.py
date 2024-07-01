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

        # Initialize variable to store which lines in line_list have a model fit to them
        self.candidate_indices = []

        # Initialize variable to store each line obj in dict form
        self.line_dicts = []

        # Load or delete spectrum data
        if not self.fresh_start:
            self.load_spectrum_data()
        else:
            self.delete_spectrum_data()

            self.spectrum.doppler_shift = self.doppler_shift_calc()

            self.flux_calc()

            self.spectrum.save_data(emission_lines)

    
    def flux_calc(self):
        '''
            Show a series of plots where the user will select if the current emission line is noise. If it is not noise, 
            then the emission line flux is calculated
            Parameters:
                        None
            Returns:
                        None
        '''
        # Iterate through each emission line
        for i, line in enumerate(self.emission_lines.line_list):
            # Get data for each line group
            group_lam = line.group_lam
            group_mask = (self.spectrum.wavelength_data > group_lam[0] - self.spectrum.line_width) & (
                self.spectrum.wavelength_data < group_lam[len(group_lam) - 1] + self.spectrum.line_width)
            group_wavelength_data = self.spectrum.wavelength_data[group_mask]
            group_flux_data = self.spectrum.flux_data[group_mask]
            
            # See if there is airglow in the lien group
            group_airglow = self.spectrum.airglow_df[(self.spectrum.airglow_df['Central Wavelength'] >= np.min(group_wavelength_data)
                                                      ) & (self.spectrum.airglow_df['Central Wavelength'] <= np.max(group_wavelength_data))]
            
            # Line information
            rest_lam = group_lam[-1] * u.AA
            line.obs_lam = self.spectrum.doppler_shift.to(u.AA,  equivalencies=u.doppler_optical(rest_lam)).value
            w0,w1 = self.wavelength_edges(group_wavelength_data)
            sumerror = (np.sum(self.spectrum.error_data[group_mask]**2 * (w1 - w0)**2))**0.5

            # Check if a model fit was successful
            if line.model_params:
                # Create a voigt profile
                voigt_profile = self.create_voigt_profile(line.model_params)

                # Calculate profile continuum and flux
                continuum = [min(voigt_profile(group_wavelength_data)) for _ in range(len(group_wavelength_data))]
                total_sumflux = np.sum((voigt_profile(group_wavelength_data))*(w1 - w0)) 

            # If the line was used for doppler calculation -> not noise 
            if line.doppler_candidate:
                self.noise_bool_list.append(False)
            else:
                # Plot basics
                legend_params, legend_strings = [], []
                sns.set_style("darkgrid")
                sns.set_theme(rc={'axes.facecolor':'#F8F5F2'})
                fig = plt.figure(figsize=(14,7))
                ax = fig.add_subplot()
                plt.title(f"Flux vs Wavelength for {line.ion}")
                fig.suptitle("Click 'y' if is noise, 'n' if not", fontweight='bold')
                plt.xlabel('Wavelength (Å)', fontsize =12)
                plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=12)
                cid = fig.canvas.mpl_connect('key_press_event', lambda event: self.on_key(event, 'Noise Detection'))
                plt.plot(group_wavelength_data, group_flux_data, linewidth=1, color = '#4B3C30')

                # Plot emission lines
                for wavelength in group_lam:
                    rest_lam = plt.axvline(x = wavelength, color = "#71816D", linewidth = 1, linestyle=((0, (5, 5))))
                obs_lam = plt.axvline(x = line.obs_lam, color = "#D7816A", linewidth = 1)

                # Check if a Voigt fit was successful
                if line.model_params:
                    # Plot voigt profile
                    voigt_fit, = plt.plot(group_wavelength_data, voigt_profile(group_wavelength_data), color = "#231651")
                    legend_params.append(voigt_fit), legend_strings.append("Voigt Fit")
                else:
                    continuum = self.emission_lines.split_create_trendline(group_wavelength_data, group_flux_data)
                    total_sumflux = np.sum(group_flux_data*(w1-w0))

                # Plot continuum
                continuum_fit, = plt.plot(group_wavelength_data, continuum, color = "#DA667B")
                legend_params.insert(0, continuum_fit), legend_strings.insert(0, "Continuum")

                # Plot airglow if applicable
                if len(group_airglow) > 0:
                    for airglow in group_airglow['Central Wavelength']:
                        airglow_lam = plt.axvline(x = airglow, color = '#4464AD', linewidth = 1)
                        legend_params.insert(0, airglow_lam), legend_strings.insert(0, "Airglow")

                legend_params.insert(0, obs_lam), legend_strings.insert(0, "Observed Wavelength")
                legend_params.insert(0, rest_lam), legend_strings.insert(0, "Rest Wavelength")
                legend = plt.legend(legend_params, legend_strings)
                legend.get_frame().set_facecolor('white')
                plt.show()
            
            # Check if noise
            if self.noise_bool_list[i]:
                flux = sumerror * (-3)
                sumerror = 0
            else:
                # Calculate final flux
                continuum_sumflux = np.sum(continuum*(w1 - w0))
                flux = total_sumflux - continuum_sumflux

            # Update parameters
            line.noise_bool = self.noise_bool_list[i]    
            line.continuum = continuum[0] # only need to store 1 value
            line.flux_error = (flux, sumerror) # SEE IF THIS FORMAT WORKS
            self.line_dicts.append(self.emission_lines.emission_line_to_dict(line))
        
        # Add data to json
        with open(self.spectrum.emission_lines_dir, "w") as json_file:
            json.dump(self.line_dicts, json_file, indent=4)


    def doppler_shift_calc(self):
        '''
            Calculates the doppler shift given a series of candidates, and writes this value to the
            according directory
            Parameters: 
                        None
            Returns:
                        doppler_shift: average doppler shift of preselected emission line groups
        '''
        # Show doppler shift plots
        self.doppler_selection_plots()

        # Calculate doppler shift
        dv = []
        # Iterate through each doppler candidate
        for i, boolean in enumerate(self.doppler_bool_list):
            group_doppler = []
            
            # If is a doppler candidate
            if boolean: 
                # Get current emissoin line group and update it's doppler component
                curr_line = self.emission_lines.line_list[self.candidate_indices[i]]
                curr_line.doppler_candidate = True

                # Calculate doppler for each emission line in that group
                for j, rest_wavelength in enumerate(curr_line.group_lam):
                    u_rest_lam = rest_wavelength * u.AA
                    u_obs_lam = curr_line.model_params[j]['x_0'] * u.AA 
                    group_doppler.append(u_obs_lam.to(u.km/u.s,  equivalencies=u.doppler_optical(u_rest_lam)))

                dv.append(sum(group_doppler)/ len(group_doppler))

        doppler_shift = sum(dv)/len(dv)

        # Store value
        with open(self.spectrum.doppler_dir, 'a') as f:
            f.write(f"{doppler_shift.value:.3f}\n")
        
        return doppler_shift

    
    def doppler_selection_plots(self):
        '''
            Shows a series of plots which the user will use to determine if an emission line is viable to be used
            for the final doppler shift calculation
            Parameters: 
                        None
            Returns:
                        None
        '''
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
                
                # Check model type -> maybe use this logic when implement Gaussian
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


    def create_voigt_profile(self, model_params):
        '''
            Creates a Voigt1D profile from a given set of parameters from a pre-fitted profile
            Parameters: 
                        model_params: Voigt1D model parameters     
            Returns:
                        composite_model: Voigt1D profile
        '''
        voigt_profiles = []

        # Iterate through all parameters
        for params in model_params:
            voigt_profile = Voigt1D(
                x_0 = params['x_0'],
                amplitude_L = params['amplitude_L'],
                fwhm_L = params['fwhm_L'],
                fwhm_G = params['fwhm_G']
            )
            voigt_profiles.append(voigt_profile)

        # Add models together
        composite_model = voigt_profiles[0]
        for voigt_profile in voigt_profiles[1:]:
            composite_model += voigt_profile

        return composite_model
    

    def wavelength_edges(self, group_wavelength_data):
        """
            Creates "boxes" to be useed to measure the flux
            Parameters: 
                        None       
            Returns:
                        w0: left edge of the boxes
                        w1: right edge of the boxes
        """
        diff = np.diff(group_wavelength_data)
        diff0 = np.concatenate((np.array([diff[0]]), diff)) 
        diff1 = np.concatenate((diff, np.array([diff[-1]]))) 
        w0 = group_wavelength_data - diff0/2.
        w1 = group_wavelength_data + diff1/2.

        return w0, w1
    

    def load_spectrum_data(self):
        '''
            Load presaved spectrum data
            Parameters: 
                        None     
            Returns:
                        None
        '''
        # Load saved doppler shift
        self.spectrum.doppler_shift = np.loadtxt(self.spectrum.doppler_dir)*(u.km/u.s)

        # Read emission_line data from JSON file
        with open(self.spectrum.emission_lines_dir, "r") as json_file:
            emission_line_data = json.load(json_file)

        # Reconstruct emission_line objects from dictionaries
        self.emission_lines.line_list.clear() 
        for data in emission_line_data:
            emission_line_obj = self.emission_lines.Emission_Line(**data) 
            self.emission_lines.line_list.append(emission_line_obj)


    def delete_spectrum_data(self):
        '''
            Delete all presaved spectrum data
            Parameters: 
                        None     
            Returns:
                        None
        '''
        try:
            os.remove(self.spectrum.doppler_dir)
            os.remove(self.spectrum.emission_lines_dir)
            os.remove(self.spectrum.fits_dir)
            os.remove(self.spectrum.ecsv_dir)
            os.remove(self.spectrum.csv_dir)
            os.remove(self.spectrum.final_plot_dir)
        except FileNotFoundError as e:
            print(e)
    

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



    