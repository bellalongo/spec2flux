# spec2flux


## About spec2flux
Spec2flux (Spectrum to Flux) aims to accurately calculate the emission line flux of lines in the Far Ultraviolet (FUV) range. This range contains information that can help clue in on exoplanetary atmospheres with the measurements providing insight into FUV radiation from host stars. Integrating UV spectra with X-ray data allows us to estimate the stellar corona, adding to the toolkit of exoplanet atmosphere data. 


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

   
2. **Create a `main.py`, which will used to put the package attributes.**

   
3. **Ensure your files for the DEM lines, airglow, and the spectrum which will used to calculate the flux are in your directory.**
   ```bash
      project-directory/
      ├── main.py
      ├── DEM_goodlinelist.csv
      ├── airglow.csv
      ├── spectrum.fits
      ```

   
4. **Make your `main.py` resemble the following:**
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

   
5. **Run `main.py`.**
   ```sh
      python main.py
   ```

   
6. **Select "best" doppler shift candidates.** </br>
   A series of plots will appear, click 'y' if the plot is an eligible candidate for doppler shift, 'n' if it is not.
   
   <div align="center">
       <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/doppler_calc.png?raw=true" alt="Doppler calculation example" width="600">
       <p>Figure 1: This is a good example because the Doppler shift is consistent between both emission lines, with minimal noise.</p>
   </div>


7. **After iterating through all the doppler shift candidates, plots will appear for the noise selection portion.** </br>
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

8. **After all emission lines are iterated through, a final plot will appear**

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

9. **After running the script, your directory will now look like the following:**
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
   
