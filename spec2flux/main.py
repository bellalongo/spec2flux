import matplotlib.pyplot as plt


from input_check import *
from emission_lines import *
from flux_calculator import *
from spectrum_data import *


def main():
    # Spectrum details (adjust me with each star)
    spectrum_dir = 'spectra/hlps'
    rest_dir = 'DEM_goodlinelist.csv'
    observation = 'sci' # SCI only
    telescope = 'hst' # HST only
    instrument = 'COS' # STIS or COS only
    grating = 'g130m' # L or M grating only
    star_name = 'toi-260'
    min_wavelength = 1160

    # Spectrum adjustments
    apply_smoothing = False # True if want to apply gaussian smoothing
    line_fit_model = 'Voigt' # 'Voigt' or 'Gaussian' fit

    # User adjustable parameters
    fresh_start = False # True if first time running, or False if have already ran for a star and want to see final plot

    # Check inputs
    InputCheck(spectrum_dir, rest_dir, observation, telescope, instrument, grating, star_name, 
               min_wavelength, apply_smoothing, line_fit_model, fresh_start)

    # Load spectrum data and emission lines
    spectrum = SpectrumData(spectrum_dir, rest_dir, observation, telescope, instrument, grating, star_name, 
                            min_wavelength, apply_smoothing)
    emission_lines = EmissionLines(spectrum)

    # Calculate flux
    flux_calc = FluxCalculator(spectrum, emission_lines, fresh_start, line_fit_model)

    # Show final plot
    spectrum.final_spectrum_plot(emission_lines, flux_calc)


if __name__ == '__main__':
    main()

    print('All done!')