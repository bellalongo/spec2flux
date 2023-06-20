from flux_calc import *
import astropy.units as u


"""
    Calculates the doppler shift based off of peaks and high liklihood rest lam lines
    Name:       doppler_shift_calc()
    Parameters: 
                rest_lam_data: rest lam data
                wavelength_data: masked wavelength data from the spectra
                flux_range: range to measure the flux of each peak
    Returns:
                doppler_shift: doppler shift of the spectra
"""
def doppler_shift_calc(rest_lam_data, wavelength_data, flux_range):
    # Initialize variables
    count = 0
    iterations = 0
    dv = []
    previous_wavelength = 0 *u.AA

    # Iterate through all of the rest wavelengths
    for wavelength in rest_lam_data["Wavelength"]:
        # Check if the current wavelength has a high liklihood to measure
        if((wavelength > 1400) and rest_lam_data['Likelihood to measure'][count] == "High"):     
            wavelength = wavelength*u.AA
            # See if there is a peak that is close to the current rest wavelength
            for peak in  wavelength_data:
                if check_in_range(peak-1.75, peak+1.75, wavelength.value):
                    # Calculate the doppler shift
                    u_rest_lam = wavelength
                    u_obs_lam = peak * u.AA
                    dv.append(u_obs_lam.to(u.km/u.s,  equivalencies=u.doppler_optical(u_rest_lam)))

                    # Check for doublet
                    blended_line_bool = blended_line_check(previous_wavelength, wavelength, iterations, flux_range)
                    if blended_line_bool:
                        # Combine the last two calculations
                        dv_length = len(dv)
                        dv[dv_length - 2] = (dv[dv_length - 2] + dv[dv_length - 1])/2
                        del dv[dv_length-1]

                    previous_wavelength = wavelength
                    break
            iterations+=1
        count+=1
    
    doppler_shift = sum(dv)/len(dv)
    return doppler_shift
