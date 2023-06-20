import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

"""
    Creates "boxes" to be useed to measure the flux
    Name:       wavelength_edges()
    Parameters: 
                wavelength_data: masked wavelength data from the spectra       
    Returns:
                w0: left edge of the box
                w1: right edge of the box
"""
def wavelength_edges(wavelength_data):
    diff = np.diff(wavelength_data)
    diff0 = np.concatenate((np.array([diff[0]]), diff)) 
    diff1 = np.concatenate((diff, np.array([diff[-1]]))) 
    w0 = wavelength_data - diff0/2.
    w1 = wavelength_data + diff1/2.

    return w0, w1


"""
    Create a continuum trendline using the average flux from the left and right of the peak
    Name:       split_create_trendline()
    Parameters: 
                wavelength_data: masked wavelength data from the spectra     
                flux_data: masked flux data from the spectra
                blended_line_bool: is the current peak blended or not
                peak_width_pixels: peak width in pixels  
    Returns:
                p(wavelength_array): continuum data for the current peak
"""
def split_create_trendline(wavelength_data, flux_data, blended_line_bool, peak_width_pixels):
    # Initialize variables
    length = len(wavelength_data) - 1
    flux_list_left = []
    flux_list_right = []
    
    # Make an array of all flux that aren't included in the peak
    for i in range(0, peak_width_pixels):
        flux_list_left.append(flux_data[i])
        flux_list_right.append(flux_data[length-i])
        
    # Find the average flux for the left and right
    avg_flux_left = sum(flux_list_left)/len(flux_list_left)
    avg_flux_right = sum(flux_list_right)/len(flux_list_right)
    
    # Use the lesser of the two values as the average flux
    if avg_flux_left < avg_flux_right:
        avg_flux = avg_flux_left
    else:
        avg_flux = avg_flux_right
    
    flux_array = [avg_flux for i in range(length + 1)]
           
    wavelength_array = np.array(wavelength_data)
    
    z = np.polyfit(wavelength_data, flux_array, 1)
    p = np.poly1d(z)
    
    return p(wavelength_array)


"""
    Create a continuum trendline using the average flux of absoption flux
    Name:       noisy_trendline()
    Parameters: 
                wavelength_data: masked wavelength data from the spectra     
                flux_data: masked flux data from the spectra
                blended_line_bool: is the current peak blended or not
    Returns:
                p(wavelength_array): continuum data for the absorption lines
"""
def noisy_trendline(wavelength_data, flux_data, blended_line_bool):  
    wavelength_array = np.array(wavelength_data)
    
    z = np.polyfit(wavelength_data, flux_data, 2)
    p = np.poly1d(z)
    
    if blended_line_bool == True:
        plt.plot(wavelength_array, p(wavelength_array), color="maroon", alpha=0.5)
    else:
        plt.plot(wavelength_array, p(wavelength_array), color="powderblue", alpha=0.5)
       
    return p(wavelength_array) 
    

"""
    Checks if an integer is within the range of two values
    Name:       check_in_range()
    Parameters: 
                val1: leftmost value   
                val2: rigthtmost value
                x: value to be compared
                peak_width_pixels: peak width in pixels  
    Returns:
                boolean: True if it is within the range, False if it is not
"""
def check_in_range(val1, val2, x):
    if val1 <= x <= val2:
        return True
    else:
        return False
    

"""
    Checks if the current peak is blended
    Name:       blended_line_check()
    Parameters: 
                previous_obs: previous observable wavlength   
                obs_lam: current observable wavelength
                iterations: number of iterations to determine the flux
                flux_range: range to measure the flux of each peak 
    Returns:
                boolean: True if the peak is blended, False if it is not
"""
def blended_line_check(previous_obs, obs_lam, iterations, flux_range):
    if check_in_range(obs_lam.value - flux_range, obs_lam.value + flux_range, previous_obs.value) and (iterations != 0):
        return True
    else:
        return False
    
