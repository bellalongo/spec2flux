from config import SPECTRUM_CONFIG, ANALYSIS_CONFIG

import sys
sys.path.insert(0, '../')

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