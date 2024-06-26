import matplotlib.pyplot as plt


from spectrum_data import *


def new_main():
    print('Loading new_main.py')

    # Spectrum details (adjust me with each star)
    spectrum_dir = 'spectra/ex_fits_file.fits'
    rest_dir = 'DEM_goodlinelist.csv'
    instrument = 'stis'
    grating = 'e140m'
    star_name = 'example'
    min_wavelength = 1160

    # User adjustable parameters
    apply_smoothing = False # True if want to apply gaussian smoothing
    fresh_start = True # True if first time running, or have already ran for a star and want to see final plot

    # Load spectrum data
    spectrum = SpectrumData(spectrum_dir, rest_dir, instrument, grating, star_name, min_wavelength, apply_smoothing)

    # next function should take in spectrum as an object

    # Sanity check
    plt.plot(spectrum.wavelength_data, spectrum.flux_data)
    plt.show()


if __name__ == '__main__':
    new_main()

    print('All done!')