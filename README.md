# Emission-Line-Flux-Calc 

### Description:
Given a fits spectrum file, the script will find the emission lines and calculate the flux for each emission line. The data is outputted as a .fits file, which can be used for analysis. 

### How to run:
Make sure the .fits file is in the same directory. <br />
Run with: **_python ./main.py 'path to fits file' 'instrument' 'grating' 'name of star'_** <br />
So for the current example, you would run with: <br />
**_python .\main.py .\spectra\hlsp_muscles_hst_stis_tau_ceti_e140m_v1_component-spec.fits 'sci' 'e140m' 'tau ceti'_** <br />

### Determining noise:
A plot will appear of each emission line before the flux is calculated, and you will click 'y' on your keyboard for noise and 'n' for not noise. </br>
![alt text](https://github.com/bellalongo/Flux-and-Line-Measurement-Script/blob/main/not_noise_example.png?raw=true)

### Using the data:
Fits Header: The header of the fits file contains all of the necessary information used to calculate the flux. <br />
Fits Data: Contains the emission line's wavelength, flux, error, and a boolean value representing if the emission line was blended. <br />
Note: Noise is marked as having the negative of the upper limit (3*error) as the flux and an error of 0. <br />

