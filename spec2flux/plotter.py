# plotting.py
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import seaborn as sns


# ------------------------------
# Plotter Class
# ------------------------------
class Plotter:
    def __init__(self, spectrum, emission_lines):
        """
            
        """
        self.spectrum = spectrum
        self.emission_lines = emission_lines
        
        # Intitialize selected lines to be empty
        self.selected_lines = {} 

        self.fig = None
        self.ax = None
        self.done_button = None

    # ------------------------------
    # Private Helper Functions
    # ------------------------------
    def _setup_figure(self, title, xlabel, ylabel, show_instructions=True):
        """

        """
        sns.set_style("whitegrid")
        self.fig, self.ax = plt.subplots(figsize=(15, 8))
        
        self.ax.set_title(title, fontsize=14)
        self.ax.set_xlabel(xlabel, fontsize=12)
        self.ax.set_ylabel(ylabel, fontsize=12)
        
        if show_instructions:
            self.fig.text(0.02, 0.98, 
                    "Click on lines to select/deselect for doppler calculation\n" +
                    "Pink = Not Selected, Red = Selected\n" +
                    "Use matplotlib toolbar to zoom/pan",
                    fontsize=10, verticalalignment='top')
        
        return self.fig, self.ax

    def _setup_interactions(self, fig):
        """

        """
        # Add Done button
        ax_done = plt.axes([0.81, 0.01, 0.1, 0.04])
        self.done_button = Button(ax_done, 'Done', color='lightgray')
        self.done_button.on_clicked(self.on_done)

        # Connect picker event
        self.fig.canvas.mpl_connect('pick_event', self._on_pick)

    def on_done(self, event):
        """

        """
        print("\nSelected lines:")
        for wavelength, data in self.selected_lines.items():
            if data['selected']:
                print(f"Line at wavelength {wavelength}: Ion = {data['ion']}")
        plt.close(self.fig)

    def _plot_spectrum(self, ax):
        """

        """
        ax.plot(self.spectrum.wavelength_data, self.spectrum.flux_data,
                linewidth=1, color='gray', alpha=0.7)
        ax.grid(True, alpha=0.3)

    def _plot_models_and_lines(self, ax, fitted_models):
        """

        """
        for group, data in fitted_models.items():
            fitted_model = data['model']
            mask = data['mask']
            ion = data['ion']

            # Plot the fitted model
            ax.plot(self.spectrum.wavelength_data[mask], 
                   fitted_model(self.spectrum.wavelength_data[mask]), 
                   color="#111D4A", alpha=0.8)
            
            # Plot the selectable lines for each wavelength
            for wavelength in group:
                # Create a wider invisible line for easier clicking
                click_line = ax.axvline(x=wavelength, color='none', linewidth=10, picker=True)
                # Create the visible line
                visible_line = ax.axvline(x=wavelength, color='pink', linestyle='--', 
                                        linewidth=2, alpha=0.8)
                
                self.selected_lines[wavelength] = {
                    'click_line': click_line,
                    'visible_line': visible_line,
                    'selected': False,
                    'ion': ion
                }

    def _interactive_plot(self, title, xlabel, ylabel, show_instructions=True):
        """

        """
        # Setup the plot
        self.fig, self.ax = self._setup_figure(title, xlabel, ylabel, show_instructions)
        
        # Plot base spectrum
        self._plot_spectrum(self.ax)
        
        # Setup interactive elements
        self._setup_interactions(self.fig)
        
        return self.fig, self.ax

    def _on_pick(self, event):
        """

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
                new_color = 'red' if new_state else 'pink'
                self.selected_lines[wavelength]['visible_line'].set_color(new_color)
        
        # Force redraw
        self.fig.canvas.draw_idle()

    def doppler_plots(self, fitted_models):
        """

        """
        # Initialize selected lines to be empty
        self.selected_lines = {}
        
        # Setup interactive plot with specific labels
        fig, ax = self._interactive_plot(
            'Emission Line Selection for Doppler Calculation',
            'Wavelength (Å)',
            'Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)'
        )
        
        # Add models and interactive lines
        self._plot_models_and_lines(ax, fitted_models)
        
        plt.show()

    def get_selected_lines(self):
        """
            
        """
        return {wavelength: self.selected_lines[wavelength]['ion'] 
                for wavelength in self.selected_lines 
                if self.selected_lines[wavelength]['selected']}
    
    
    # def _plot_base(self, title, xlabel, ylabel, figtext):
    #     """

    #     """
    #     sns.set_style("whitegrid")
    #     fig = plt.figure(figsize=(14, 7))
    #     ax = fig.add_subplot()

    #     plt.title(title, fontsize = 14)
    #     plt.xlabel(xlabel, fontsize = 12)
    #     plt.ylabel(ylabel, fontsize = 12)
        
    #     if figtext:
    #         fig.text(0.02, 0.98, 
    #                  "Click on lines to select/deselect for doppler calculation\n" +
    #                  "Pink = Not Selected, Red = Selected\n" +
    #                  "Use matplotlib toolbar to zoom/pan",
    #                  fontsize=10, verticalalignment='top')

    #     return fig, ax

    # def _add_legend(self, ax, handles, labels):
    #     """
            
    #     """
    #     legend = ax.legend(handles, labels)
    #     legend.get_frame().set_facecolor('white')
    #     return legend
    
    # def _interactive_plot(self, title, xlabel, ylabel):
    #     """

    #     """
    #     fig, ax = self._plot_base(self, title, xlabel, ylabel)

    #     # Plot the full spectrum
    #     plt.plot(self.spectrum.wavelength_data, self.spectrum.flux_data,
    #              linewidth = 1, color = 'gray', alpha = 0.7)
        
    #     # Add Done button
    #     ax_done = plt.axes([0.85, 0.02, 0.1, 0.04])
    #     btn_done = Button(ax_done, 'Done')
    #     btn_done.on_clicked(lambda event: plt.close(fig))

    #     # Connect picker event
    #     fig.canvas.mpl_connect('pick_event', self._on_pick)

    #     # Enable zoom toolbar
    #     plt.grid(True, alpha=0.3)

    # def _on_pick(self, event):
    #     """

    #     """
    #     curr_line = event.artist

    #     # Find which wavelength was clicked and get its group
    #     clicked_wavelength = None
    #     clicked_group = None

    #     # Find the clicked wavelengths
    #     for wavelength in self.selcted_lines:
    #         if self.selected_lines[wavelength]['click_line'] == curr_line:
    #             clicked_wavelength = wavelength
    #             break

    #     if clicked_wavelength is None:
    #         return
        
    #     # Find which group this wavelegth belongs to
    #     for ion in self.line_list:
    #         for emission_line in self.line_list[ion].values():
    #             if clicked_wavelength in emission_line.group_lam:
    #                 clicked_group = emission_line.group_lam
    #                 break
        
    #     if clicked_group is None:
    #         return
        
    #     # Toggle selection state based on clicked line
    #     new_state = not self.selected_lines[clicked_wavelength]['selected']

    #     # Update all lines in the group
    #     for wavelength in clicked_group:
    #         self.selected_lines[wavelength]['selected'] = new_state

    #         # Update color of visible line
    #         new_color = 'red' if new_state else 'pink' # UPDATE ME
    #         self.selected_lines[wavelength]['visible_line'].set_color(new_color)

    #     # Force redraw
    #     curr_line.figure.canvas.draw_idle()

    # def _get_selected_lines(self):
    #     """

    #     """
    #     return {
    #         wavelength: self.selected_lines[wavelength] 
    #         for wavelength in self.selected_lines
    #         if self.selected_lines[wavelength]['selected']
    #     }
 
    # # ------------------------------
    # # Public Methods
    # # ------------------------------
    # def doppler_plots(self, fitted_models):
    #     """

    #     """
    #     # Initialize selected lines to be empty
    #     self.selected_lines = {} 

    #     # Setup plot bones -> maybe fix me
    #     self._interactive_plot(
    #         'Emission Line Selection for Doppler Calculation',
    #         'Wavelength (Å)',
    #         'Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)') 
        
    #     # Iterate through fitted models
    #     for group, data in fitted_models.items():
    #         fitted_model = data['model']
    #         mask = data['mask']
    #         ion = data['ion']

    #         # Plot the fitted model
    #         plt.plot(self.spectrum.wavelength_data[mask], 
    #                 fitted_model(self.spectrum.wavelength_data[mask]), 
    #                 color="#111D4A", alpha=0.8)
            
    #         # Plot the selectable lines for each wavelength
    #         for wavelength in group:
    #             # Create a wider invisible line for easier clicking
    #             click_line = plt.axvline(x=wavelength, color='none', linewidth=10, picker=True)

    #             # Create the visible line
    #             visible_line = plt.axvline(x=wavelength, color='pink', linestyle='--', 
    #                                      linewidth=2, alpha=0.8)
                
    #             self.selected_lines[wavelength] = {
    #                 'click_line': click_line,
    #                 'visible_line': visible_line,
    #                 'selected': False,
    #                 'ion': ion
    #             }

    #     plt.show()
        
        
        



        
        # # Plot fitted models and lines for each group
        # for group_data in fitted_models.items():

        #     ax.plot(self.spectrum.wavelength_data[data['mask']], 
        #            data['model'](self.spectrum.wavelength_data[data['mask']]), 
        #            color="#111D4A", alpha=0.8)

        


    # def noise_plot(self, line, wavelength_data, flux_data, group_airglow, model_profile=None, continuum=None):
    #     """

    #     """
    #     # Create base plot
    #     fig, ax = self._plot_base(
    #         title=f"Flux vs Wavelength for {line.ion}",
    #         xlabel='Wavelength (Å)',
    #         ylabel='Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)'
    #     )
    #     fig.suptitle("Click 'y' if is noise, 'n' if not", fontweight='bold', fontsize=15)
        
    #     # Plot data
    #     ax.plot(wavelength_data, flux_data, linewidth=1, color='#4B3C30')
    #     legend_handles = []
    #     legend_labels = []

    #     # Add rest/observed wavelengths
    #     rest_lam = ax.axvline(x=line.group_lam[-1], color="#71816D", linewidth=1, linestyle=((0, (5, 5))))
    #     obs_lam = ax.axvline(x=line.obs_lam, color="#D7816A", linewidth=1)
    #     legend_handles.extend([rest_lam, obs_lam])
    #     legend_labels.extend(["Rest Wavelength", "Observed Wavelength"])

    #     # Add model fit
    #     if model_profile is not None:
    #         model_fit, = ax.plot(wavelength_data, model_profile(wavelength_data), color="#231651")
    #         legend_handles.append(model_fit)
    #         legend_labels.append(f"Model Fit")

    #     # Add continuum
    #     if continuum is not None:
    #         continuum_fit, = ax.plot(wavelength_data, continuum, color="#DA667B")
    #         legend_handles.append(continuum_fit)
    #         legend_labels.append("Continuum")

    #     # Add airglow
    #     if not group_airglow.empty:
    #         for airglow in group_airglow['Central Wavelength']:
    #             airglow_lam = ax.axvline(x=airglow, color='#4464AD', linewidth=1)
    #         legend_handles.append(airglow_lam)
    #         legend_labels.append("Airglow")

    #     self._add_legend(ax, legend_handles, legend_labels)
    #     return fig

    # def doppler_selection_plot(self, ion, group_mask, fitted_model, group):
    #     """

    #     """
    #     # Create base plot
    #     fig, ax = self._plot_base(
    #         title=f"Flux vs Wavelength for {ion}",
    #         xlabel='Wavelength (Å)',
    #         ylabel='Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)'
    #     )
    #     fig.suptitle("Click 'y' to select line for Doppler calculation", fontweight='bold', fontsize=15)

    #     # Plot data and model
    #     ax.plot(self.spectrum.wavelength_data[group_mask], self.spectrum.flux_data[group_mask], linewidth=1)
    #     ax.plot(self.spectrum.wavelength_data[group_mask], fitted_model(self.spectrum.wavelength_data[group_mask]), color="#111D4A")

    #     # Add interactive lines
    #     legend_handles = []
    #     legend_labels = []
    #     for wavelength in group:
    #         rest_lam = ax.axvline(x=wavelength, color="#F96E46", linestyle=((0, (5, 5))), linewidth=1)
    #         obs_lam = ax.axvline(x=wavelength, color="#8F1423", linewidth=1)
    #     legend_handles.extend([rest_lam, obs_lam])
    #     legend_labels.extend(["Rest Wavelength", "Observed Wavelength"])

    #     self._add_legend(ax, legend_handles, legend_labels)
    #     return fig