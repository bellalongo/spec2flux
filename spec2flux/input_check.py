from os.path import exists


class InputCheck(object):
    def __init__(self, spectrum_config, analysis_config):

        # Extract from configs
        self.spectrum_dir = spectrum_config['spectrum_dir']
        self.rest_dir = spectrum_config['rest_dir']
        self.airglow_dir = spectrum_config['airglow_dir']
        self.observation = spectrum_config['observation'].upper()
        self.telescope = spectrum_config['telescope'].upper()
        self.instrument = spectrum_config['instrument'].upper()
        self.grating = spectrum_config['grating'].upper()
        self.resolution = spectrum_config['resolution'].upper()
        self.star_name = spectrum_config['star_name']
        self.min_wavelength = spectrum_config['min_wavelength']
        self.max_wavelength = spectrum_config['max_wavelength']

        self.apply_smoothing = analysis_config['apply_smoothing']
        self.line_fit_model = analysis_config['line_fit_model']
        self.cont_fit = analysis_config['cont_fit']
        self.fresh_start = analysis_config['fresh_start']

        # Check files
        self.check_files()

        # Check spectrum data
        self.check_spectrum_data()

        # Check booleans
        self.check_booleans()

        # Check continuum fit
        self.check_continuum()

        # Check model
        self.check_model()


    def check_files(self):
        """
            Checks if files exist and if are of correct variable type
            Parameters: 
                        None
            Returns:
                        None
        """
        # Check types
        if not isinstance(self.spectrum_dir, str):
            raise TypeError(f"Variable spectrum_dir must be of type 'str'")
        
        if not isinstance(self.rest_dir, str):
            raise TypeError(f"Variable rest_dir must be of type 'str'")
        
        if not isinstance(self.airglow_dir, str):
            raise TypeError(f"Variable airglow_dir must be of type 'str'")

        # Check values
        if not exists(self.spectrum_dir):
            raise FileNotFoundError(f"The file for spectrum_dir: {self.spectrum_dir} does not exist")
        
        if not exists(self.rest_dir):
            raise FileNotFoundError(f"The file for rest_dir: {self.rest_dir} does not exist")
        
        if not exists(self.airglow_dir):
            raise FileNotFoundError(f"The file for rest_dir: {self.airglow_dir} does not exist")
        
    
    def check_spectrum_data(self):
        """
            Checks if spectrum data inputs are of correct variable types and value
            Parameters: 
                        None
            Returns:
                        None
        """
        # Check types
        if not isinstance(self.observation, str):
            raise TypeError(f"Variable observation must be of type 'str'")       
        
        if not isinstance(self.telescope, str):
            raise TypeError(f"Variable telescope must be of type 'str'")

        if not isinstance(self.instrument, str):
            raise TypeError(f"Variable instrument must be of type 'str'")
        
        if not isinstance(self.grating, str):
            raise TypeError(f"Variable grating must be of type 'str'")
        
        if not isinstance(self.resolution, str):
            raise TypeError(f"Variable resolution must be of type 'str'")
        
        if not isinstance(self.star_name, str):
            raise TypeError(f"Variable star_name must be of type 'str'")
        
        if not isinstance(self.min_wavelength, (int, float)):
            print(self.min_wavelength, type(self.min_wavelength))
            raise TypeError(f"Variable min_wavelength must be of type 'int' or 'float'")
        
        if not isinstance(self.max_wavelength, (int, float)):
            print(self.max_wavelength, type(self.max_wavelength))
            raise TypeError(f"Variable min_wavelength must be of type 'int' or 'float'")
        
        # Check values
        if self.observation != 'SCI':
            raise ValueError(f"Observation {self.observation} is not compatible")
        
        if self.telescope != 'HST':
            raise ValueError(f"Telescope {self.telescope} is not compatible")

        if self.instrument != 'STIS' and self.instrument != 'COS':
            raise ValueError(f"Instrument {self.instrument} is not compatible")
        
        if 'L' not in self.grating and 'M' not in self.grating:
            raise ValueError(f"Grating {self.grating} is not compatible")
        
        if self.resolution != 'HIGH' and self.resolution != 'LOW':
            raise ValueError(f"Resolution {self.resolution} is not compatible")
        
    
    def check_booleans(self):
        """
            Checks if input booleans are of type bool
            Parameters: 
                        None
            Returns:
                        None
        """
        # Check types
        if not isinstance(self.apply_smoothing, bool):
            raise TypeError(f"Variable apply_smoothing must be of type 'bool'")
        
        if not isinstance(self.fresh_start, bool):
            raise TypeError(f"Variable fresh_start must be of type 'bool'")
        

    def check_continuum(self):
        """
            Checks if the specified continuum fit is of correct type and value
            Parameters:
                        None
            Returns:
                        None
        """
        # Check types
        if not isinstance(self.cont_fit, str):
            raise TypeError(f"Variable cont_fit must be of type 'str'")
        
        # Check value
        if self.cont_fit != 'Complete' and self.cont_fit != 'Individual':
            raise ValueError(f"Continuum fit {self.cont_fit} is not compatible")


    def check_model(self):
        """
            Checks if model fit variable is of correct type and value
            Parameters: 
                        None
            Returns:
                        None
        """
        # Check type
        if not isinstance(self.line_fit_model, str):
            raise TypeError(f"Variable line_fit_model must be of type 'str'")

        # Check value
        if self.line_fit_model != 'Voigt' and self.line_fit_model != 'Gaussian':
            raise ValueError(f"Model type {self.line_fit_model} is not compatible")