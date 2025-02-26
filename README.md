# spec2flux

[![astropy](https://img.shields.io/badge/astropy-red?style=for-the-badge)](https://astropy.org/)
[![matplotlib](https://img.shields.io/badge/matplotlib-orange?style=for-the-badge)](https://matplotlib.org/stable/index.html)
[![scipy](https://img.shields.io/badge/scipy-yellow?style=for-the-badge)](https://scipy.org/)
[![seaborn](https://img.shields.io/badge/seaborn-green?style=for-the-badge)](https://seaborn.pydata.org/)
[![numpy](https://img.shields.io/badge/numpy-blue?style=for-the-badge)](https://numpy.org/doc/)
[![pandas](https://img.shields.io/badge/pandas-purple?style=for-the-badge)](https://pandas.pydata.org/docs/)

## Table of Contents
- [About](#about-spec2flux)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Interactive Analysis Workflow](#interactive-analysis-workflow)
- [Output Files](#output-files)

## About spec2flux

Spec2flux (Spectrum to Flux) accurately calculates emission line flux in the Far Ultraviolet (FUV) range. These measurements provide critical insights into FUV radiation from exoplanet host stars, helping researchers understand exoplanetary atmospheres. When integrated with X-ray data, these measurements allow estimation of stellar corona properties.

### How It Works

The script performs flux calculations through the following steps:

1. **Emission Line Grouping**: Groups emission lines using a configurable tolerance
2. **Doppler Shift Calculation**: Fits composite Voigt profiles to strong emission lines and calculates Doppler shift
3. **Flux Calculation**: Processes each line, filtering noise, and calculating flux using either Voigt profiles or raw data
4. **Results Generation**: Creates publication-ready plots and exports calculations as FITS and ECSV files

## Installation

```bash
pip install spec2flux
```

## Configuration

Spec2flux uses two configuration dictionaries to control its behavior:

### SPECTRUM_CONFIG

```python
SPECTRUM_CONFIG = {
    'spectrum_dir': 'spectrum-file.fits',       # Path to spectrum FITS file
    'rest_dir': 'DEM_goodlinelist.csv',         # Path to emission line REST file
    'airglow_dir': 'airglow.csv',               # Path to airglow CSV file
    'observation': 'sci',                       # Observation type ('sci' only)
    'telescope': 'hst',                         # Telescope used ('hst' only)
    'instrument': 'stis',                       # Instrument ('stis' or 'cos')
    'grating': 'e140m',                         # Grating type (must contain 'L' or 'M')
    'resolution': 'high',                       # Resolution ('high' or 'low')
    'star_name': 'STAR_NAME',                   # Target star name (used for output files)
    'min_wavelength': 1160,                     # Minimum wavelength for analysis
    'max_wavelength': 1700                      # Maximum wavelength for analysis
}
```

### ANALYSIS_CONFIG

```python
ANALYSIS_CONFIG = {
    'apply_smoothing': False,                   # Apply gaussian smoothing
    'line_fit_model': 'Voigt',                  # Fitting model ('Voigt' or 'Gaussian')
    'cont_fit': 'Individual',                   # Continuum fit method ('Complete' or 'Individual')
    'fresh_start': True                         # True for first run or to regenerate plots
}
```

## Usage

Create a `main.py` file that uses the spec2flux toolkit:

```python
import spec2flux
from config import SPECTRUM_CONFIG, ANALYSIS_CONFIG

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
        flux_calculator = spec2flux.FluxCalculator(spectrum, emission_lines, plotter)
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

Create a `config.py` file with your configuration settings (see [Configuration](#configuration) section).

Your project directory should be organized as follows:

```
project-directory/
├── main.py
├── config.py
├── DEM_goodlinelist.csv
├── airglow.csv
├── spectrum.fits
```

Run your analysis:

```bash
python main.py
```

## Interactive Analysis Workflow

<details>
<summary><b>Click to expand interactive workflow details</b></summary>

### 1. Doppler Shift Candidate Selection

When the script runs, you'll first be presented with a spectrum plot for selecting Doppler shift candidates:

<p align="center">
  <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/doppler_spectrum.png?raw=true" alt="Doppler selection" width="700">
</p>

**Steps:**
1. Use the zoom button to inspect potential candidates (orange lines)
2. Click on good candidates to select them (they will turn magenta) **zoom button must be unclicked to select!**
3. Use the home button to return to the full spectrum view
4. Click "Done" when you've selected all suitable Doppler candidates

<p align="center">
  <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/doppler_candidates_selected.png?raw=true" alt="Completed doppler selection" width="700">
</p>

### 2. Emission Line Selection

After Doppler selection, you'll see the spectrum with potential emission line candidates:

<p align="center">
  <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/noise_spectrum.png?raw=true" alt="Emission line selection" width="700">
</p>

**Steps:**
1. Use the zoom button to inspect orange emission line candidates
2. Click on real emission lines to select them (they will turn magenta)
3. Leave noise as orange
4. Click "Done" when you've identified all emission lines

### 3. Final Results

After completing the selection process, a final plot will be generated showing the full spectrum with the continuum, profile fits, and marked wavelengths:

<p align="center">
  <img src="https://github.com/bellalongo/spec2flux/blob/main/readme_pics/final_plot.png?raw=true" alt="Final analysis plot" width="700">
</p>

</details>

## Output Files

After running the analysis, your directory will contain new files organized as follows:

```
project-directory/
├── ... (original files)
├── doppler/
│   ├── star_name_doppler.txt           # Doppler shift calculation results
├── emission_lines/
│   ├── star_name_lines.json            # Detailed emission line data
├── flux/
│   ├── star_name.csv                   # Flux calculations (CSV format)
│   ├── star_name.ecsv                  # Flux calculations (ECSV format)
│   ├── star_name.fits                  # Flux calculations (FITS format)
├── plots/
│   ├── star_name_final_plot.png        # Final spectrum plot with all components
```

### Output Data Format

The output files contain the following information:

#### Header Information
* `DATE`: Date flux was calculated
* `FILENAME`: Name of the FITS file used for flux calculation
* `TELESCP`: Telescope used to measure spectrum
* `INSTRMNT`: Active instrument used to measure spectrum
* `GRATING`: Grating used to measure spectrum
* `TARGNAME`: Name of target star
* `DOPPLER`: Doppler shift used to measure flux
* `WIDTH`: Average peak width of emission lines
* `RANGE`: Flux range used to isolate emission line
* `WIDTHPXL`: Average emission line peak width in pixels
* `UPPRLIMIT`: Upper limit used for noise (3*error)

#### Data Columns
* `Ion`: Ion identifier
* `Rest Wavelength`: Rest wavelength of the emission line
* `Flux`: Calculated flux value
* `Error`: Error estimate for the flux
* `Blended Line`: Boolean flag indicating blended lines

**Note**: Noise is marked with the negative of the upper limit (3*error) as the flux and an error of 0.
