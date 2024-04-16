import numpy as np
import astropy.units as u
from scipy.ndimage import gaussian_filter

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
    Smooths data using gaussian_filter
    Name:       smooth_data()
    Parameters: 
                wavelength_data: masked wavelength data from the spectra  
                flux_data: masked flux data from the spectra   
                error_data: masked error data from the spectra  
                sigma: standard deviation of the gaussian (controls how much neighboring data points contribute to smoothing)
    Returns:
                smoothed_wavelength: smoothed wavelength data
                smoothed_flux: smoothed flux data
                smoothed_error: smoothed error data
"""
def smooth_data(wavelength_data, flux_data, error_data, sigma):
    # Smooth data
    smoothed_wavelength = gaussian_filter(wavelength_data, sigma)
    smoothed_flux = gaussian_filter(flux_data, sigma)
    smoothed_error = gaussian_filter(error_data, sigma)

    return smoothed_wavelength, smoothed_flux, smoothed_error


"""
    Create a continuum trendline using the average flux from the left and right of the peak
    Name:       split_create_trendline()
    Parameters: 
                wavelength_data: masked wavelength data from the spectra     
                flux_data: masked flux data from the spectra
                peak_width_pixels: peak width in pixels  
    Returns:
                p(wavelength_array): continuum data for the current peak
"""
def split_create_trendline(wavelength_data, flux_data, peak_width_pixels):
    # Initialize variables
    length = len(wavelength_data) - 1
    flux_list_left = []
    flux_list_right = []
    
    # Make an array of all flux that aren't included in the peak
    for i in range(0, int(peak_width_pixels/2)):
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
    
    continuum_array = [avg_flux for i in range(length + 1)]
    
    return continuum_array
    
