import numpy as np

def subchannel_splice(waves,specs,dspecs):
    """
    Splice together spectra from multiple MIRI subchannels.
    """

    #Sort out our arrays based on wavelength so that we can extract the scaling factors
    wave_srt, spec_srt, dspec_srt = zip(*sorted(zip(waves,specs,dspecs), key = lambda x: x[0]))

    wave_shrt  = np.array(wave_srt[0])
    spec_shrt  = np.array(spec_srt[0])
    dspec_shrt = np.array(dspec_srt[0])

    wave_med  = np.array(wave_srt[1])
    spec_med  = np.array(spec_srt[1])
    dspec_med = np.array(dspec_srt[1])

    wave_lng  = np.array(wave_srt[2])
    spec_lng  = np.array(spec_srt[2])
    dspec_lng = np.array(dspec_srt[2])

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



    # #Find the overlap between SHORT and MEDIUM channels
    # shrt_ind = np.where(wave_shrt >= min(wave_med))
    # med_shrt_ind = np.where(wave_med <= max(wave_shrt))

    # #Find the median of SHORT and MEDIUM in the overlapping ranges
    # shrt_median = np.nanmedian(spec_shrt[shrt_ind])
    # med_shrt_median = np.nanmedian(spec_med[med_shrt_ind])
    # shrt_factor = med_shrt_median / shrt_median

    # #Find the SHORT region without overlap
    # shrt_exc = np.where(wave_shrt <= min(wave_med))

    # #Scale the SHORT spectrum and concatenate to the MEDIUM spectrum
    # shrt_spec_scaled = spec_shrt * shrt_factor
    # shrt_dspec_scaled = dspec_shrt * shrt_factor
    
    # wave_shrt_med = np.concatenate([wave_shrt[shrt_exc],wave_med])
    # spec_shrt_med = np.concatenate([shrt_spec_scaled[shrt_exc],spec_med])
    # dspec_shrt_med = np.concatenate([shrt_dspec_scaled[shrt_exc],dspec_med])

    # #Find the overlap between the MEDIUM and LONG sub-channels
    # lng_ind = np.where(wave_lng <= max(wave_shrt_med))
    # med_lng_ind = np.where(wave_shrt_med >= min(wave_lng))

    # #Find the median of MEDIUM and LONG in the overlapping ranges
    # lng_median     = np.nanmedian(spec_lng[lng_ind])
    # med_lng_median = np.nanmedian(spec_shrt_med[med_lng_ind])
    # lng_factor = med_lng_median / lng_median

    # #Find the LNG region without overlap
    # lng_exc = np.where(wave_lng >= max(wave_shrt_med))

    # #Scale the LNG spectrum and concatenate it to the SHT-MED spectrum
    # lng_spec_scaled  = spec_lng * lng_factor
    # lng_dspec_scaled = dspec_lng * lng_factor

    # wave_shrt_med_lng    = np.concatenate([wave_shrt_med,wave_lng[lng_exc]])
    # spec_shrt_med_lng  = np.concatenate([spec_shrt_med,lng_spec_scaled[lng_exc]])
    # dspec_shrt_med_lng = np.concatenate([dspec_shrt_med,lng_dspec_scaled[lng_exc]])

    # return wave_shrt_med_lng, spec_shrt_med_lng, dspec_shrt_med_lng