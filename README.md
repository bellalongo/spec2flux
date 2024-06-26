# fluxpal


## About fluxpal
The goal of fluxpal is to accurately calculate the emission line flux of lines in the Far Ultraviolet(FUV) range. This range contains information that can help clue in on exoplanetary atmospheres with the measurements providing insight into FUV radiation from host stars. Integrating UV spectra with X-ray data allows us to estimate the stellar corona, adding to the toolkit of exoplanet atmosphere data. 


## How it works
The script accomplishes the task of calculating flux by:
1. Grouping emission lines using a preset tolerance, which can be adjusted if using a different DEM (CSV) file.
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


## Installation
1. Clone the repository </br>

   ```sh
   git clone https://github.com/bellalongo/fluxpal.git
   ```
2. Move the fits file being used to calculate flux into the 'spectra' folder.
3. Adjust the 'fresh_start' parameter depending on whether you want to adjust calculations on a new run </br>
   (if first run, set to True)
   ```sh
   fresh_start = True
   ```
3. Adjust the 'gaussian_smoothing' parameter depending on whether you to smooth the spectrum </br>
   (if first run, set to False)
   ```sh
   gaussian_smoothing = False
   ```
5. Run the script:
    ```sh
   python main.py spectra/ex_fits_file.fits 'instrument' 'grating' 'star name'
   ```

## Adjustments
The script can be run on its own without any adjustments, but if the star is a little finicky:
### main.py
* adjust the mask used to isolate emission lines
* edit header if want different header information
### emission_lines.py
* adjust peak_width in 'peak_width_finder'
* adjust tolerance used to group emission lines in 'grouping_emission_lines'
* adjust voigt fit initial parameters
  


## Doppler shift calculation
A plot will appear of the 'best' emission lines, with a Voigt fit fitted to it. Click 'y' if you think the line should be used to calculate Doppler shift (these lines will automatically not be considered as noise), and 'n' if the emission line if not. </br>

![doppler calculation example](https://github.com/bellalongo/Flux-and-Line-Measurement-Script/blob/main/readme_pics/doppler_calc.png?raw=true)
</br>

## Determining if noise
All lines not selected for Doppler calculation will appear, some with a Voigt profile fitted if possible to be selected as noise or not noise. Click 'y' if you think the line is noise, and 'n' if the line is not noise. </br>

**Not noise:**
![Not noise example](https://github.com/bellalongo/Flux-and-Line-Measurement-Script/blob/main/readme_pics/not_noise.png?raw=true)
</br>
</br>
**Noise:**
![Not noise example](https://github.com/bellalongo/Flux-and-Line-Measurement-Script/blob/main/readme_pics/noise.png?raw=true)

## Final plot
After all emission line selections and calculations have been made, a final plot will appear showing the spectrum and each labeled emission line. Matplotlib gives the ability to zoom into the plot, so please do so to double-check the lines. This plot is saved to the 'plots' folder. </br>

![final plot example](https://github.com/bellalongo/Flux-and-Line-Measurement-Script/blob/main/plots/example_final_plot.png?raw=true)
</br>



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


### Using the data:
ECSV Data: File columns are Ion, Rest Wavelength, Flux, Error, Blended Line, each row representing a grouped emission line.  <br />
FITS Data: FITS Table containing the same data as the ECSV file. <br />
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
