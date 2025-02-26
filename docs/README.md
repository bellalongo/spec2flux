# spec2flux - Developer Guide

This README provides instructions for running the code directly from the cloned repository.

## Getting Started

1. Clone the repository
   ```sh
   git clone https://github.com/bellalongo/spec2flux.git
   ```

2. Navigate to the docs directory
   ```sh
   cd spec2flux/docs
   ```

3. Ensure your spectrum file is available in the docs directory

4. Open `tutorial.py` and adjust configuration settings:
   - Update the SPECTRUM_CONFIG dictionary with your spectrum details
   - Set analysis parameters in ANALYSIS_CONFIG

5. Set the appropriate run parameters:
   ```python
   # Set to True for first run or if you want to regenerate all outputs
   # Set to False to use previously saved data
   'fresh_start': True
   
   # Set to True if you want to apply Gaussian smoothing to the spectrum
   # Generally recommended to leave as False for first analysis
   'apply_smoothing': False
   ```

6. Run the script:
   ```sh
   python tutorial.py
   ```

## Adjustments for Troubleshooting

If you're having issues getting good line measurements, you can adjust several parameters in the code:

### In tutorial.py
* **Wavelength range**: Adjust `min_wavelength` and `max_wavelength` values if your emission lines are being cut off
  ```python
  'min_wavelength': 1160,  # Minimum wavelength to analyze
  'max_wavelength': 1700   # Maximum wavelength to analyze
  ```

* **Line fit model**: Switch between 'Voigt' and 'Gaussian' to see which works better for your data
  ```python
  'line_fit_model': 'Voigt',  # Try 'Gaussian' if Voigt profiles aren't fitting well
  ```

* **Continuum fitting**: Change between 'Complete' and 'Individual' methods
  ```python
  'cont_fit': 'Individual',  # Try 'Complete' for challenging spectra
  ```

### In emission_lines.py
* **Tolerance**: Adjust the grouping tolerance for emission lines
  ```python
  # Near the top of the file:
  TOLERANCE = 6  # Tolerance for grouping emission lines (in Angstroms)
  ```

* **Peak width finder**: In the `peak_width_finder` method, you can adjust the width used to isolate emission lines
  ```python
  # Increase for broader lines, decrease for narrower lines
  peak_width = PEAK_WIDTH_LOW_RES if self.resolution == 'LOW' else PEAK_WIDTH_HIGH_RES
  ```

* **Voigt fit parameters**: Adjust initial parameters for better fits
  ```python
  # In _create_voigt_profile method:
  init_amp = np.max(self.spectrum.flux_data[line_mask])
  init_fwhm_g = self.spectrum.line_width/5  # Try adjusting this divisor
  init_fwhm_l = self.spectrum.line_width/5  # Try adjusting this divisor
  ```

### In spectrum_data.py
* **Peak widths**: Modify the constants for high and low resolution spectra
  ```python
  # Near the top of the file:
  PEAK_WIDTH_LOW_RES = 5.0  # Angstroms
  PEAK_WIDTH_HIGH_RES = 0.5  # Angstroms
  ```

* **Smoothing**: Adjust the smoothing sigma if using gaussian smoothing
  ```python
  # In smooth_data method, adjust the sigma parameter:
  def smooth_data(self, sigma):
      # Default is sigma=1, try values between 0.5-2 for different smoothing levels
  ```

### In flux_calculator.py
* **Upper limit factor**: Change the factor used for calculating upper limits from errors
  ```python
  # Near the top of the file:
  UPPER_LIMIT_FACTOR = 3  # Factor for calculating upper limit from error
  ```

## Output Directory Structure

After running the analysis, the following directories will be created:

```
docs/
├── doppler/                         # Doppler shift calculations
│   └── star_name_doppler.txt
├── emission_lines/                  # Emission line data
│   └── star_name_lines.json
├── flux/                            # Flux calculation results
│   ├── star_name.csv
│   ├── star_name.ecsv
│   └── star_name.fits
└── plots/                           # Generated plots
    └── star_name_final_plot.png
```
