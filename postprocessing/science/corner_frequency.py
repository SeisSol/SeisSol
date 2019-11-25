#!/usr/bin/env python3

#Python (Version 3.6) functions for determining the Spectral Seismological misfit, corner frequency and falloff rate 
#of a seismological amplitude spectra
#
#by Nico Schliwa

import numpy as np
#T. Oliphant. A guide to numpy, 2006. URL http://www.numpy.org/

#Function for linear interpolation

def linear_interpolate(amplitude_trace, old_time_comp, start_frequency, end_frequency, frequency_steps):
    """Usually spectra do not have equal frequency steps, but this is needed for further functions of this script.
    With this function a spectrum gets linear interpolated and the time compiler (frequency steps) 
    transformed to equal density."""
    
    interpolate_range=np.arange(start_frequency, end_frequency+frequency_steps, frequency_steps)
    new_amplitude_trace=np.zeros(len(interpolate_range))

    for q in range(0,len(interpolate_range)):
        for d in range(0,len(old_time_comp)):
            if old_time_comp[d] > interpolate_range[q]:
                new_amplitude_trace[q]=abs(amplitude_trace[d-1])+(interpolate_range[q]-old_time_comp[d-1])*(
                                        abs(amplitude_trace[d])-abs(amplitude_trace[d-1]))/(
                                        old_time_comp[d]-old_time_comp[d-1])
                break

    new_time_comp = interpolate_range
    
    return new_amplitude_trace, new_time_comp


def spectral_seismological_misfit(first_amplitude_trace, second_amplitude_trace):
    """The spectral seismological misfit quantifies the misfit between two spectra.
    The time compilers (frequency steps) need to be identical. The differences to other misfits
    are that single outliers are not weighted as strong as for example by the RMSD and relative
    misfits are not scaled by their amplitude. Definition after:
    S. Karimzadeh, A. Askan, and A. Yakut. Assessment of Simulated Ground Motions in 
    Earthquake Engineering Practice: A Case Study for Duzce (Turkey), 2017"""
    
    misfit = np.mean(abs(np.log10(abs(first_amplitude_trace)/abs(second_amplitude_trace))))
    
    return misfit


def spectral_corner_frequency(amplitude_trace, time_comp, startfreq=0., 
                             endfreq=0., misfit="SSM", n=2, plot=False):
    """Function to determine the corner frequency of a seismic amplitude spectrum after the Brune model 
    for shear waves in the far-field (J. N. Brune. Tectonic stress and the spectra of seismic shear waves 
    from earthquakes, 1970). It calculates the function A/(1+(f/fc)^n) for every corner frequency fc in the
    chosen interval and returns the fc, for which the misfit between the function and the spectrum is 
    minimal. The amplitude A is the amplitude of the second lowest frequency of the spectrum (for more details:
    B. P. Allmann, P. M. Shearer. Global variations of stress drop for moderate to large earthquakes, 2009).
    In the classic Brune model the falloff rate is n=2, but there are several other models with different 
    falloff rates, therefore the falloff rate n can be manually changed. Start and endpoint define the frequency 
    interval, in which the corner frequency is searched. Default is the whole time compiler, but reducing the 
    interval can make the calculation much faster. The Spectral Seismological misfit is used by default, with
    misfit="RMSD" it can be changed to the root-mean-square deviation. By turning the plot parameter to TRUE 
    as second returning the amplitudes array of the whole plot, which was used to minimize the misfit is given."""
    
    if endfreq==0.: 
        endfreq=time_comp[len(time_comp)-1]
    if startfreq==0.:
        startfreq=time_comp[1]
    amplitude=amplitude_trace[1]
    corner_function=np.zeros_like(amplitude_trace)
    corner_bestfit=9999999.0
    
    for fc in np.arange(startfreq, endfreq, time_comp[1]-time_comp[0]):
        
        for x in range(0,len(time_comp),1):
            corner_function[x]=amplitude*(1/(1+(time_comp[x]/fc)**n))
        if misfit=="SSM":
            corner_misfit=np.mean(abs(np.log10(abs(corner_function)/abs(amplitude_trace))))
        if misfit=="RMSD":
            corner_misfit=np.sqrt(np.mean((corner_function-amplitude_trace)**2))
        if corner_misfit<corner_bestfit:
                corner_bestfit=corner_misfit
                bestfit_corner_frequency=fc
                bestfit_corner_function=amplitude*(1/(1+(time_comp/fc)**n))
    
    if plot:
        return bestfit_corner_frequency, bestfit_corner_function
    else:
        return bestfit_corner_frequency

    
    
def spectral_corner_frequency_and_falloff(amplitude_trace, time_comp, startfreq=0., 
                             endfreq=0., startrate=1.5, endrate=3.5, misfit="SSM", plot=False):
    """Nearly the same function as spectral_corner_frequency, but here the falloff rate n is also
    a variable that is searched for to minimize the misfit. Considered n are between the start- and 
    endrate in steps of 0.1. This function is computationally very expensive."""
    
    if endfreq==0.: 
        endfreq=time_comp[len(time_comp)-1]
    if startfreq==0.:
        startfreq=time_comp[1]
    amplitude=amplitude_trace[1]
    corner_function=np.zeros_like(amplitude_trace)
    corner_bestfit=9999999.0
    
    for fr in np.arange(startrate, endrate, 0.1):
        for fc in np.arange(startfreq, endfreq, time_comp[1]-time_comp[0]):

            for x in range(0,len(time_comp),1):
                corner_function[x]=amplitude*(1/(1+(time_comp[x]/fc)**fr))
            if misfit=="SSM":
                corner_misfit=np.mean(abs(np.log10(abs(corner_function)/abs(amplitude_trace))))
            if misfit=="RMSD":
                corner_misfit=np.sqrt(np.mean((corner_function-amplitude_trace)**2))
            if corner_misfit<corner_bestfit:
                    corner_bestfit=corner_misfit
                    bestfit_corner_frequency=fc
                    bestfit_falloff_rate=fr
                    bestfit_corner_function=amplitude*(1/(1+(time_comp/fc)**fr))

    if plot:
        return bestfit_corner_frequency, bestfit_falloff_rate, bestfit_corner_function
    else:
        return bestfit_corner_frequency, bestfit_falloff_rate
