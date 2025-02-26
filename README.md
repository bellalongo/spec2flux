# spec2flux


## About spec2flux
Spec2flux (Spectrum to Flux) aims to accurately calculate the emission line flux of lines in the Far Ultraviolet (FUV) range. This range contains information that can help clue in on exoplanetary atmospheres, with the measurements providing insight into FUV radiation from host stars. Integrating UV spectra with X-ray data allows us to estimate the stellar corona, adding to the toolkit of exoplanet atmosphere data. 


## How it works
The script accomplishes the task of calculating flux by:
1. Grouping emission lines using a preset tolerance, which can be adjusted if using a different DEM file.
2. Calculating the Doppler shift by fitting a composite Voigt profile to the strongest emission lines, and comparing Voigt peaks to rest wavelength to calculate.
3. Iterating through each line, determining if noise, and then calculating the flux either using the Voigt profile or the original data.
4. Presents a final plot of all of the final emission lines, and exports the calculations as a FITS and ECSV file.


## Prerequisites 
* [![astropy][astropy-pic]][astropy-url]
* [![matplotlib][matplotlib-pic]][matplotlib-url]
* [![scipy][scipy-pic]][scipy-url]
* [![seaborn][seaborn-pic]][seaborn-url]
* [![numpy][numpy-pic]][numpy-url]
* [![pandas][pandas-pic]][pandas-url]

</br>
</br>


## Pip Installation
1. **Pip install the repository in your desired directory.**
   
   ```sh
      pip install spec2flux
   ```
   
2. **Create a `main.py`, which will be utilized to perform the calculations using the spec2flux toolkit.**
   ```python
      import spec2flux
            
      def main():
          # Check inputs
          spec2flux.InputCheck(SPECTRUM_CONFIG, ANALYSIS_CONFIG)
      
          # Load spectrum data and emission lines
          spectrum = spec2flux.SpectrumData(SPECTRUM_CONFIG, ANALYSIS_CONFIG)
          emission_lines = spec2flux.EmissionLines(spectrum)
      
          if ANALYSIS_CONFIG['fresh_start']:
              # Initialize plotter
              plotter = spec2flux.Plotter(spectrum, emission_lines)
      
              # Calculate doppler shift
              plotter.doppler_plots()
              emission_lines.update_selected_lines(plotter, 'Doppler')
      
              doppler_calculator = spec2flux.DopplerCalculator(spectrum)
              spectrum.doppler_shift = doppler_calculator.calculate_doppler_shift(emission_lines)
      
              # Handle noise detection
              plotter.noise_plots()
              emission_lines.update_selected_lines(plotter, 'Noise')
      
              # Calculate fluxes
              flux_calculator = spec2flux.FluxCalculator(spectrum, emission_lines)
              flux_calculator.process_spectrum()
      
          else:
              # Load existing data
              emission_lines.load_saved_data()
              
          # Final plot
          plotter = spec2flux.Plotter(spectrum, emission_lines)
          plotter.create_final_plot()
      
      if __name__ == '__main__':
          main()
   ```

3. **Create a `config.py`, which contains two main configuration dictionaries that control the behavior of spec2flux.**
   ``` python
   """
       Configuration settings for spectrum details and preferences
   
           spectrum_dir: Path to spectrum FITS file
           rest_dir: Path to emission line REST file
           airglow_dir: Path to airglow CSV file
           observation: Observation type ('sci' only)
           telescope: Telescope used ('hst' only)
           instrument: Instrument used ('stis' or 'cos')
           grating: Grating type
           resolution: High or low resolution ('high' or 'low')
           star_name: Name of target star
           min/max_wavelength: Wavelength range to analyze
      """
      SPECTRUM_CONFIG = {
          'spectrum_dir': 'hlsp_muscles_hst_stis_tau_ceti_e140m_v1_component-spec.fits',
          'rest_dir': 'DEM_goodlinelist.csv',
          'airglow_dir': 'airglow.csv',
          'observation': 'sci',
          'telescope': 'hst',
          'instrument': 'stis',
          'grating': 'e140m',
          'resolution': 'high',
          'star_name': 'TESTING',
          'min_wavelength': 1160,
          'max_wavelength': 1700
      }
      
      """
          Configuration settings for emission line fitting and analysis    
      
              apply_smoothing: Apply gaussian smoothing
              line_fit_model: Fitting model ('Voigt' or 'Gaussian')
              cont_fit: Continuum fitting method ('Complete' or 'Individual')
              fresh_start: True for first run or to regenerate plots
      """
      ANALYSIS_CONFIG = {
          'apply_smoothing': False,
          'line_fit_model': 'Voigt',
          'cont_fit': 'Individual',
          'fresh_start': True
      }
   ```

4. **Ensure your files for the DEM lines, airglow, and the spectrum, which will used to calculate the flux, are in your directory.**
   ```bash
      project-directory/
      ├── main.py
      ├── config.py
      ├── DEM_goodlinelist.csv
      ├── airglow.csv
      ├── spectrum.fits
      ```

5. **Make your `main.py` resemble the following:**
   ```python
      import spec2flux
      
      
      def main():
          # Spectrum details (adjust me with each star)
          spectrum_dir = 'spectrum-spec.fits'
          rest_dir = 'DEM_goodlinelist.csv'
          airglow_dir = 'airglow.csv'
          observation = 'sci' # SCI only
          telescope = 'hst' # HST only
          instrument = 'stis' # STIS or COS only
          grating = 'e140m' # L or M grating only
          star_name = 'NEWEX'
          min_wavelength = 1160
      
          # Spectrum adjustments
          apply_smoothing = False # True if want to apply gaussian smoothing
          line_fit_model = 'Voigt' # 'Voigt' or 'Gaussian' fit
      
          # User adjustable parameters
          fresh_start = True # True if first time running, or have already ran for a star and want to see final plot
      
          # Check inputs
          spec2flux.InputCheck(spectrum_dir, rest_dir, airglow_dir, observation, telescope, instrument, grating, star_name, 
                     min_wavelength, apply_smoothing, line_fit_model, fresh_start)
      
          # Load spectrum data and emission lines
          spectrum = spec2flux.SpectrumData(spectrum_dir, rest_dir, airglow_dir, observation, telescope, instrument, grating, star_name, 
                                  min_wavelength, apply_smoothing)
          emission_lines = spec2flux.EmissionLines(spectrum)
      
          # Calculate flux
          flux_calc = spec2flux.FluxCalculator(spectrum, emission_lines, fresh_start, line_fit_model)
      
          # Show final plot
          spectrum.final_spectrum_plot(emission_lines, flux_calc)
      
      
      if __name__ == '__main__':
          main()
   ```

   
6. **Run `main.py`.**
   ```sh
      python main.py
   ```

7. **Select "best" doppler shift candidates.** </br>
A plot of the complete spectrum will appear, with orange lines representing possible candidates for doppler shift.
 

   
7. **Select "best" doppler shift candidates.** </br>
   A series of plots will appear, click 'y' if the plot is an eligible candidate for doppler shift, 'n' if it is not.
   
   <div align="center">
       <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/doppler_calc.png?raw=true" alt="Doppler calculation example" width="600">
       <p>Figure 1: This is a good example because the Doppler shift is consistent between both emission lines, with minimal noise.</p>
   </div>


8. **After iterating through all the doppler shift candidates, plots will appear for the noise selection portion.** </br>
   To select if the current plot is noise or not, click 'y' if the plot is noise and 'n' if the plot is not noise.

   <table>
     <tr>
       <td align="center">
         <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/not_noise.png?raw=true" alt="Not noise example" width="400">
         <p>Figure 2: This plot is not noisy because the emission lines are well-defined at the specified rest wavelengths.</p>
       </td>
       <td align="center">
         <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/noise.png?raw=true" alt="Noise example" width="400">
         <p>Figure 3: This plot is noisy because the emission lines at the specified rest wavelength are not well-defined.</p>
       </td>
     </tr>
   </table>

9. **After all emission lines are iterated through, a final plot will appear**

   <div align="center">
       <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/final_plot.png?raw=true" alt="Final plot example" width="600">
       <p>Figure 4: Plot of the entire spectrum with the continuum, profile fits, and wavelengths marked.</p>
   </div>
   
   </br>
   <div align="center">
       <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/zoom.png?raw=true" alt="Zoomed plot example" width="600">
       <p>Figure 5: Final plot zoomed in.</p>
   </div>
   </br>

10. **After running the script, your directory will now look like the following:**
      ```bash
      project-directory/
      ├── main.py
      ├── DEM_goodlinelist.csv
      ├── airglow.csv
      ├── spectrum.fits
      ├── doppler/
      │   ├── newex_doppler.txt
      ├── emission_lines/
      │   ├── newex_lines.json
      ├── flux/
      │   ├── newex.csv
          ├── newex.ecsv
          ├── newex.fits
      ├── plots/
      │   ├── newex_final_plot.png
      ```

## Using the data:
Header contains:
* DATE: date flux was calculated
* FILENAME: name of the fits file used to for flux calc
* TELESCP: telescope used to measure spectrum
* INSTRMNT: active instrument to measure spectrum
* GRATING: grating used to measure spectrum
* TARGNAME: name of star used in measurement
* DOPPLER: doppler shift used to measure flux
* WIDTH: average peak width of the emissoin lines
* RANGE: flux range used to isolate emission line
* WIDTHPXL: average emission line peak width in pixels
* UPPRLIMIT: upper limit used for noise

ECSV + CSV Data: File columns are Ion, Rest Wavelength, Flux, Error, Blended Line, each row representing a grouped emission line.  <br />
FITS Data: FITS Table contains the same data as the ECSV file. <br />
Note: Noise is marked as having the negative of the upper limit (3*error) as the flux and an error of 0. <br />



[astropy-url]: https://astropy.org/
[astropy-pic]: https://img.shields.io/badge/astropy-red?style=for-the-badge
[matplotlib-url]: https://matplotlib.org/stable/index.html
[matplotlib-pic]: https://img.shields.io/badge/matplotlib-orange?style=for-the-badge
[scipy-url]: https://scipy.org/
[scipy-pic]: https://img.shields.io/badge/scipy-yellow?style=for-the-badge
[seaborn-url]: https://seaborn.pydata.org/
[seaborn-pic]: https://img.shields.io/badge/seaborn-green?style=for-the-badge
[numpy-url]: https://numpy.org/doc/
[numpy-pic]: https://img.shields.io/badge/numpy-blue?style=for-the-badge
[pandas-url]: https://pandas.pydata.org/docs/
[pandas-pic]: https://img.shields.io/badge/pandas-purple?style=for-the-badge
   
