import numpy as np

def subchannel_splice(waves,specs,dspecs):
    """
    Splice together spectra from multiple MIRI subchannels. Sorted by ascending wavelength of each segment.

    Inputs
        waves - list containing the wavelengths for each segment to be spliced together
        specs - list containing the spectra for each segment to be spliced together
        dspecs - list containing the uncertainties for each segment to be spliced together

    Outputs
        lists for the spliced wavelengths, spectra, and uncertainties

    """

    #Sort out our arrays based on wavelength so that we can extract the scaling factors
    wave_srt, spec_srt, dspec_srt = zip(*sorted(zip(waves,specs,dspecs), key = lambda x: x[0]))

    wave = []
    spec = []
    dspec = []

    for i in range(len(waves)):
        if i==0:
            wave = np.concatenate([wave,np.array(wave_srt[0])])
            spec = np.concatenate([spec,np.array(spec_srt[0])])
            dspec = np.concatenate([dspec,np.array(dspec_srt[0])])
        else:
            #Find the region of overlap between the concatenated spectrum and the next region
            low_ind = np.where(wave >= min(wave_srt[i]))
            high_ind = np.where(wave_srt[i] <= max(wave))

            #Find the median in the overlapping regions
            low_median = np.nanmedian(np.array(spec)[low_ind])
            high_median = np.nanmedian(np.array(spec_srt[i])[high_ind])
            factor = low_median / high_median

            #Find the region of the concatenated spectrum without overlap
            exc = np.where(wave <= min(wave_srt[i]))

            #Scale the new region and concatenate
            scale_spec = factor * np.array(spec_srt[i])
            scale_dspec = factor * np.array(dspec_srt[i])

            wave = np.concatenate([wave[exc],wave_srt[i]])
            spec = np.concatenate([spec[exc],scale_spec])
            dspec = np.concatenate([dspec[exc],scale_dspec])

    return wave, spec, dspec


