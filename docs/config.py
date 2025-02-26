"""
    Configuration settings for spectrum details and preferences

        spectrum_dir: Path to spectrum FITS file
        rest_dir: Path to emission line REST file
        airglow_dir: Path to airglow CSV file
        observation: Observation type ('sci' only)
        telescope: Telescope used ('hst' only)
        instrument: Instrument used ('stis' or 'cos')
        grating: Grating type
        resolution: High or low resolution ('high' or 'low')
        star_name: Name of target star
        min/max_wavelength: Wavelength range to analyze
"""
SPECTRUM_CONFIG = {
    'spectrum_dir': 'hlsp_muscles_hst_stis_tau_ceti_e140m_v1_component-spec.fits',
    'rest_dir': 'DEM_goodlinelist.csv',
    'airglow_dir': 'airglow.csv',
    'observation': 'sci',
    'telescope': 'hst',
    'instrument': 'stis',
    'grating': 'e140m',
    'resolution': 'high',
    'star_name': 'TESTING',
    'min_wavelength': 1160,
    'max_wavelength': 1700
}

"""
    Configuration settings for emission line fitting and analysis    

        apply_smoothing: Apply gaussian smoothing
        line_fit_model: Fitting model ('Voigt' or 'Gaussian')
        cont_fit: Continuum fitting method ('Complete' or 'Individual')
        fresh_start: True for first run or to regenerate plots
"""
ANALYSIS_CONFIG = {
    'apply_smoothing': False,
    'line_fit_model': 'Voigt',
    'cont_fit': 'Individual',
    'fresh_start': True
}