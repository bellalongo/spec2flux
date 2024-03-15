# fluxpal


## About fluxpal
The goal of fluxpal is to accurately calculate the emission line flux of lines in the Far Ultraviolet(FUV) range. This range contains information that can help clue in on exoplanetary atmospheres with the measurements providing insight into FUV radiation from host stars. Integrating UV spectra with X-ray data allows us to estimate the stellar corona, adding to the toolkit of exoplanet atmosphere data. 


## How it works
The script accomplishes the task of calculating flux by:
1. Grouping emission lines using a preset tolerance, which can be adjusted if using a different DEM (CSV) file.
2. Calculating the doppler shift by fitting a composite Voigt profile to the strongest emission lines, and comparing Voigt peaks to rest wavelength to calculate.
3. Iterating through each line, determing if noise, and then calculating the flux either using the Voigt profile or the original data.
4. Presents a final plot of all of the final emission lines, and exports the calculations as a FITS and ECSV file.


## Packages
* [![astropy][astropy-pic]][astropy-url]
* [![matplotlib][matplotlib-pic]][matplotlib-url]
* [![seaborn][seaborn-pic]][seaborn-url]
* [![numpy][numpy-pic]][numpy-url]
* [![pandas][pandas-pic]][pandas-url]


  







[numpy-url]: https://numpy.org/doc/
[numpy-pic]: https://img.shields.io/badge/numpy-blue?style=for-the-badge
[astropy-url]: https://astropy.org/
[astropy-pic]: https://img.shields.io/badge/astropy-orange?style=for-the-badge
[matplotlib-url]: https://matplotlib.org/stable/index.html
[matplotlib-pic]: https://img.shields.io/badge/matplotlib-yellow?style=for-the-badge
[seaborn-url]: https://seaborn.pydata.org/
[seaborn-pic]: https://img.shields.io/badge/seaborn-purple?style=for-the-badge
[pandas-url]: https://pandas.pydata.org/docs/
[pandas-pic]: https://img.shields.io/badge/pandas-black?style=for-the-badge
