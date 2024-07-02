import matplotlib.pyplot as plt


from emission_line import *
from flux_calculator import *
from spectrum_data import *


def new_main():
    print('Loading new_main.py')

    # Spectrum details (adjust me with each star)
    spectrum_dir = 'spectra/hlsp_muscles_hst_stis_tau_ceti_e140m_v1_component-spec.fits'
    rest_dir = 'DEM_goodlinelist.csv'
    instrument = 'stis'
    grating = 'e140m'
    star_name = 'NEWEX'
    min_wavelength = 1160

    # Spectrum adjustments
    apply_smoothing = False # True if want to apply gaussian smoothing
    line_fit_model = 'Gaussian' # 'Voigt' or 'Gaussian' fit

    # User adjustable parameters
    fresh_start = True # True if first time running, or have already ran for a star and want to see final plot

    # Load spectrum data and emission lines
    spectrum = SpectrumData(spectrum_dir, rest_dir, instrument, grating, star_name, min_wavelength, apply_smoothing)
    emission_lines = EmissionLines(spectrum)

    # Calculate flux
    flux_calc = FluxCalculator(spectrum, emission_lines, fresh_start, line_fit_model)

    # Show final plot
    spectrum.final_spectrum_plot(emission_lines, flux_calc)


if __name__ == '__main__':
    new_main()

    print('All done!')