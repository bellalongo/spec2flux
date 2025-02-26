from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import seaborn as sns

from .model_fitter import ModelFitter


# ------------------------------
# Plotter Class
# ------------------------------
class Plotter:
    """
        This class handles visualization of spectra and emission lines, 
        including interactive plotting for user selection of emission lines for
        Doppler calculations and noise detection.
        Attributes:
            spectrum (SpectrumData): The spectrum data object containing wavelength and flux data
            emission_lines (EmissionLines): Object containing emission line data
            model_fitter (ModelFitter): Object for fitting models to emission lines

            selected_lines (dict): Dictionary tracking user-selected emission lines
            fig (Figure): Matplotlib figure object for the current plot
            ax (Axes): Matplotlib axes object for the current plot
            done_button (Button): Button widget for completing the selection process
    """
    def __init__(self, spectrum, emission_lines):
        """
            Initializes the Plotter object with spectrum and emission line data.
            Arguments:
                spectrum (SpectrumData): The spectrum data object
                emission_lines (EmissionLines): Object containing emission line data       
        """
        self.spectrum = spectrum
        self.emission_lines = emission_lines

        self.model_fitter = ModelFitter(spectrum)
        
        # Intitialize selected lines to be empty
        self.selected_lines = {} 

        self.fig = None
        self.ax = None
        self.done_button = None

    # ------------------------------
    # Private Helper Functions
    # ------------------------------
    def _setup_figure(self, title, xlabel, ylabel, instructions):
        """
            Sets up the figure and axes for plotting with appropriate styling and labels.
            Arguments:
                title (str): Title for the plot
                xlabel (str): Label for the x-axis
                ylabel (str): Label for the y-axis
                instructions (str): Instructions to display above the plot, or None
            Returns:
                tuple: (fig, ax) Matplotlib figure and axes objects
        """
        sns.set_style("whitegrid")
        self.fig, self.ax = plt.subplots(figsize=(15, 8))
        
        self.ax.set_title(title, fontsize=16, pad = 10, fontweight = 'bold', y = -0.14)
        self.ax.set_xlabel(xlabel, fontsize=12)
        self.ax.set_ylabel(ylabel, fontsize=12)

        # Calculate plot bounds
        wavelength_mask = (self.spectrum.wavelength_data > 1250) & (
            self.spectrum.wavelength_data < self.spectrum.max_wavelength)
        max_peak = max(self.spectrum.flux_data[wavelength_mask])
        min_peak = min(self.spectrum.flux_data)

        # Set axis limits
        self.ax.set_ylim(min_peak, max_peak * 1.5)
        self.ax.set_xlim(self.spectrum.min_wavelength, self.spectrum.max_wavelength)
        
        # Add instructions as subtitle
        if instructions: 
            subtitle = (
                f"Click on lines to select/deselect for {instructions.lower()} calculation\n"
                f"Orange = {'Not Selected' if instructions == 'Doppler' else 'Noise'}, "
                f"Magenta = {'Selected' if instructions == 'Doppler' else 'Not Noise'}\n"
                "Use matplotlib toolbar to zoom/pan"
            )
                
            self.fig.suptitle(subtitle, fontsize=12, y=0.97)

        return self.fig, self.ax

    def _setup_interactions(self, fig):
        """
            Sets up interactive elements for the plot, including the Done button
            and event handlers.
            Arguments:
                fig (Figure): Matplotlib figure to add interactions to
            Returns:
                None
        """
        # Add Done button
        ax_done = plt.axes([0.81, 0.01, 0.1, 0.04])
        self.done_button = Button(ax_done, 'Done', color='lightgray')
        self.done_button.on_clicked(self._on_done)

        # Connect picker event
        self.fig.canvas.mpl_connect('pick_event', self._on_pick)

    def _on_done(self, event):
        """
            Callback function for the Done button. Prints selected lines and closes the figure.
            Arguments:
                event: The button click event
            Returns:
                None
        """
        print("\nSelected lines:")
        for wavelength, data in self.selected_lines.items():
            if data['selected']:
                print(f"Line at wavelength {wavelength}: Ion = {data['ion']}")
        plt.close(self.fig)

    def _plot_spectrum(self, ax):
        """
            Plots the basic spectrum data on the given axes.
            Arguments:
                ax (Axes): Matplotlib axes to plot on
            Returns:
                None
        """
        ax.plot(self.spectrum.wavelength_data, self.spectrum.flux_data,
                linewidth=1, color='gray', alpha=0.7)
        ax.grid(True, alpha=0.3)

    def _plot_doppler_data(self, ax):
        """
            Plots emission lines for Doppler shift calculation with interactive elements.
            Arguments:
                ax (Axes): Matplotlib axes to plot on
            Returns:
                None
        """
        # Iterate through each ion in the line list
        for ion in self.emission_lines.line_list:
            # Iterate through each emission line in this ion's group
            for emission_line in self.emission_lines.line_list[ion].values():
                # Create mask for the emission line group
                group_mask = (
                    (self.spectrum.wavelength_data > emission_line.group_lam[0] - self.spectrum.line_width) &
                    (self.spectrum.wavelength_data < emission_line.group_lam[-1] + self.spectrum.line_width)
                )
                
                # If there are model parameters, create and plot the model
                if emission_line.model_params:
                    # Create the model using the saved parameters
                    model = self.model_fitter.create_model_profile(emission_line.model_params)
                    
                    # Plot the fitted model
                    fitted_model, = ax.plot(
                        self.spectrum.wavelength_data[group_mask],
                        model(self.spectrum.wavelength_data[group_mask]),
                        color="#242325",
                        alpha=0.8
                    )
                    
                    # Plot the lines for each wavelength in the group
                    for i, wavelength in enumerate(emission_line.group_lam):
                        # Plot observed wavelength
                        if self.spectrum.line_fit_model == 'Voigt':
                            obs_wavelength = emission_line.model_params[i]['x_0']
                        else:
                            obs_wavelength = emission_line.model_params[i]['mean']
                            
                        # Calculate y-range for the observed wavelength line
                        local_flux = self.spectrum.flux_data[group_mask]
                        local_max = max(local_flux)
                        local_min = min(local_flux)
                        y_range = local_max - local_min
                        
                        # Draw line segment instead of full vertical line
                        obs_line = ax.plot(
                            [obs_wavelength, obs_wavelength],
                            [local_min - 0.1 * y_range, local_max + 0.1 * y_range],
                            color="#BBB891",
                            linewidth=1.5,
                            alpha=0.8,
                            solid_capstyle='butt'
                        )[0]

                        # Create a wider invisible line for easier clicking
                        click_line = ax.axvline(
                            x=wavelength,
                            color='none',
                            linewidth=10,
                            picker=True
                        )
                        
                        # Create the visible line
                        visible_line = ax.axvline(
                            x=wavelength,
                            color='#DC965A',
                            linestyle='--',
                            linewidth=2,
                            alpha=0.8
                        )
                        
                        # Store the line references and metadata
                        self.selected_lines[wavelength] = {
                            'click_line': click_line,
                            'visible_line': visible_line,
                            'selected': False,
                            'ion': ion,
                            'emission_line': emission_line
                        }
        proxy_visible_line = Line2D([0], [0], color='#8B1E3F', linestyle='--', linewidth=2, alpha=0.8)
        ax.legend([fitted_model, visible_line, proxy_visible_line, obs_line], ['Model Fit', 
                                                                               'Rest Wavelength (Not Doppler Candidate)', 
                                                                               'Rest Wavelength (Doppler Candidate)',
                                                                               'Observed Wavelength'])

    def _plot_noise_data(self, ax):
        """
            Plots emission lines for noise detection with interactive elements.
            Arguments:
                ax (Axes): Matplotlib axes to plot on
            Returns:
                None
        """
        # Iterate through each ion in the line list
        for ion in self.emission_lines.line_list:
            # Iterate through each emission line in this ion's group
            for emission_line in self.emission_lines.line_list[ion].values():
                # Create mask for the emission line group
                group_mask = (
                    (self.spectrum.wavelength_data > emission_line.group_lam[0] - self.spectrum.line_width) &
                    (self.spectrum.wavelength_data < emission_line.group_lam[-1] + self.spectrum.line_width)
                )

                # Calculate and plot continuum
                wavelength_data = self.spectrum.wavelength_data[group_mask]
                flux_data = self.spectrum.flux_data[group_mask]
                
                if self.spectrum.cont_fit == 'Complete':
                    continuum = [self.spectrum.spectrum_continuum] * len(wavelength_data)
                else:
                    if emission_line.model_params:
                        model_profile = self.model_fitter.create_model_profile(emission_line.model_params)
                        continuum = [min(model_profile(wavelength_data))] * len(wavelength_data)
                    else:
                        continuum = self.emission_lines._create_trendline(wavelength_data, flux_data)

                # Make sure not a doppler candidate
                if not emission_line.doppler_candidate:
                    # Plot continuum
                    continuum, = ax.plot(wavelength_data, continuum, color="#416788", lw=2, alpha=0.8)

                    # If there are model parameters, create and plot the model
                    if emission_line.model_params:
                        # Create the model using the saved parameters
                        model = self.model_fitter.create_model_profile(emission_line.model_params)
                        
                        # Plot the fitted model
                        fitted_model, = ax.plot(
                            self.spectrum.wavelength_data[group_mask],
                            model(self.spectrum.wavelength_data[group_mask]),
                            color="#242325",
                            alpha=0.8,
                            label=f"{ion} Model Fit"
                        )
                        
                    # Plot the lines for each wavelength in the group
                    for wavelength in emission_line.group_lam:
                        # Create a wider invisible line for easier clicking
                        click_line = ax.axvline(
                            x=wavelength,
                            color='none',
                            linewidth=10,
                            picker=True
                        )
                        
                        # Create the visible line
                        visible_line = ax.axvline(
                            x=wavelength,
                            color='#DC965A',
                            linestyle='--',
                            linewidth=2,
                            alpha=0.8
                        )
                        
                        # Store the line references and metadata
                        self.selected_lines[wavelength] = {
                            'click_line': click_line,
                            'visible_line': visible_line,
                            'selected': False,
                            'ion': ion,
                            'emission_line': emission_line
                        }

        proxy_visible_line = Line2D([0], [0], color='#8B1E3F', linestyle='--', linewidth=2, alpha=0.8)
        ax.legend([fitted_model, visible_line, proxy_visible_line, continuum], ['Model Fit', 
                                                                                'Rest Wavelength (Noise Candidate)',
                                                                                'Rest Wavelength (Emission Line Candidate)',
                                                                                'Continuum'])

    def _interactive_plot(self, title, xlabel, ylabel, instructions):
        """
            Creates an interactive plot with appropriate labels and setup.
            Arguments:
                title (str): Title for the plot
                xlabel (str): Label for the x-axis
                ylabel (str): Label for the y-axis
                instructions (str): Instructions to display above the plot
            Returns:
                tuple: (fig, ax) Matplotlib figure and axes objects
        """
        # Setup the plot
        self.fig, self.ax = self._setup_figure(title, xlabel, ylabel, instructions)
        
        # Plot base spectrum
        self._plot_spectrum(self.ax)
        
        # Setup interactive elements
        self._setup_interactions(self.fig)
        
        return self.fig, self.ax

    def _on_pick(self, event):
        """
            Callback function for line picking events. Toggles selection state of emission lines.
            Arguments:
                event: The picker event
            Returns:
                None
        """
        thisline = event.artist
        
        # Find which wavelength was clicked
        clicked_wavelength = None
        clicked_group = None
        
        for wavelength in self.selected_lines:
            if self.selected_lines[wavelength]['click_line'] == thisline:
                clicked_wavelength = wavelength
                break
                
        if clicked_wavelength is None:
            return
            
        # Find the group this wavelength belongs to
        for group in self.emission_lines.grouped_lines.values():
            for subgroup in group:
                if clicked_wavelength in subgroup:
                    clicked_group = subgroup
                    break
            if clicked_group:
                break
                
        if clicked_group is None:
            return
            
        # Toggle selection state
        new_state = not self.selected_lines[clicked_wavelength]['selected']
        
        # Update all lines in the group
        for wavelength in clicked_group:
            if wavelength in self.selected_lines:
                self.selected_lines[wavelength]['selected'] = new_state
                new_color = '#8B1E3F' if new_state else '#DC965A'
                self.selected_lines[wavelength]['visible_line'].set_color(new_color)
        
        # Force redraw
        self.fig.canvas.draw_idle()

    def _plot_emission_line_components(self, ax, emission_line):
        """
            Plots individual components of an emission line for the final plot.
            Arguments:
                ax (Axes): Matplotlib axes to plot on
                emission_line (EmissionLine): Emission line object to plot
            Returns:
                list: List of legend elements
        """
        legend_items = []
        
        # Create mask for the line group
        group = emission_line.group_lam
        flux_mask = (self.spectrum.wavelength_data > group[0] - self.spectrum.line_width) & (
            self.spectrum.wavelength_data < group[-1] + self.spectrum.line_width)
        
        # Plot continuum
        continuum = [emission_line.continuum] * len(self.spectrum.wavelength_data[flux_mask])
        trendline = ax.plot(
            self.spectrum.wavelength_data[flux_mask], 
            continuum, 
            color="#C05746", 
            lw=2
        )[0]
        legend_items.append({'element': trendline, 'label': 'Continuum'})

        # Plot model fit if available
        if emission_line.model_params:
            model_profile = self.model_fitter.create_model_profile(emission_line.model_params)
            model_fit = ax.plot(
                self.spectrum.wavelength_data[flux_mask],
                model_profile(self.spectrum.wavelength_data[flux_mask]),
                color="#49506F"
            )[0]
            legend_items.append({'element': model_fit, 'label': 'Model Profile'})

        # Plot rest wavelengths
        for wavelength in emission_line.group_lam:
            if emission_line.noise_bool:
                line = ax.axvline(
                    x=wavelength,
                    color='#92B257',
                    linewidth=1.5,
                    linestyle='--'
                )
                legend_items.append({'element': line, 'label': 'Noise Wavelength'})
            else:
                line = ax.axvline(
                    x=wavelength,
                    color='#5D7A3E',
                    linewidth=1.5,
                    linestyle='--'
                )
                legend_items.append({'element': line, 'label': 'Rest Wavelength'})

        return legend_items

    # ------------------------------
    # Public Methods
    # ------------------------------
    def doppler_plots(self):
        """
            Creates an interactive plot for selecting emission lines for Doppler shift calculation.
            Arguments:
                None
            Returns:
                None
        """
        # Initialize selected lines to be empty
        self.selected_lines = {}
        
        # Setup interactive plot with specific labels
        fig, ax = self._interactive_plot(
            f'{self.spectrum.star_name} Emission Line Selection for Doppler Calculation',
            'Wavelength (Å)',
            'Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)',
            'Doppler'
        )
        
        # Add models and interactive lines
        self._plot_doppler_data(ax)
        plt.show()

    def noise_plots(self):
        """
            Creates an interactive plot for selecting which emission lines are noise.
            Arguments:
                None
            Returns:
                None
        """
        # Initialize selected lines to be empty
        self.selected_lines = {}
        
        # Setup interactive plot with specific labels
        fig, ax = self._interactive_plot(
            f'{self.spectrum.star_name} Emission Line Selection for Noise Detection',
            'Wavelength (Å)',
            'Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)',
            'Noise'
        )
        
        # Add models and interactive lines
        self._plot_noise_data(ax)
        
        plt.show()

    def get_selected_lines(self):
        """
            Returns a dictionary of selected emission lines.
            Arguments:
                None
            Returns:
                dict: Dictionary of selected lines mapped to their ion names
        """
        return {wavelength: self.selected_lines[wavelength]['ion'] 
                for wavelength in self.selected_lines 
                if self.selected_lines[wavelength]['selected']}
    
    def create_final_plot(self):
        """
        Creates the final plot showing all emission lines with their models, continuum,
            and wavelength markers. Saves the plot to a file and displays it.
        Arguments:
            None
        Returns:
            None
        """
        # Setup the plot basics
        _, ax = self._setup_figure(
            title=f"Flux vs Wavelength for {self.spectrum.star_name}",
            xlabel='Wavelength (Å)',
            ylabel='Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)',
            instructions=None  # No instructions needed for final plot
        )

        # Plot base spectrum
        self._plot_spectrum(ax)

        # Initialize legend elements
        legend_elements = []
        legend_labels = []
        
        # Plot emission lines and their components
        for ion in self.emission_lines.line_list:
            for emission_line in self.emission_lines.line_list[ion].values():
                legend_items = self._plot_emission_line_components(ax, emission_line)
                
                # Add new legend items if they don't exist yet
                for item in legend_items:
                    if item['label'] not in legend_labels:
                        legend_elements.append(item['element'])
                        legend_labels.append(item['label'])

        # Add legend and save
        ax.legend(legend_elements, legend_labels)
        plt.savefig(self.spectrum.final_plot_dir)
        plt.show()