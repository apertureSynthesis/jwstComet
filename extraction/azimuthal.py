import numpy as np
import astropy.units as u
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats
from jwstComet.utils import readCube, subchannel_splice, readHeader

class azimuthal(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @u.quantity_input(waveLo = u.um, waveUp = u.um, innerRadius = u.arcsec, outerRadius = u.arcsec)
    def extractSpec(self,cubeFile,waveLo,waveUp,innerRadius=0.0*u.arcsec,outerRadius=0.1*u.arcsec,withPlots=False):
        """
        Extract a spectrum within a specified annulus, take the azimuthal average, and save it to a text file.
        Plot the aperture and extracted spectrum if desired.
        """
        try:
            if waveLo.value >= waveUp.value:
                raise ValueError('Lower wavelength must be smaller than upper wavelength')
        except (ValueError):
            exit('Update wavelength range')

        try:
            if innerRadius.value >= outerRadius.value:
                raise ValueError('Inner annulus radius must be smaller than outer annulus radius')
        except (ValueError):
            exit('Update inner/outer annulus radii')


        combined_waves = []
        combined_specs = []
        combined_sigmas = []

        for j in range(len(cubeFile)):

            #Read in the file and get relevant information and data
            sciCube = readCube(cubeFile[j])
            data = sciCube.data
            derr = sciCube.derr
            dnpts = data.shape[0]

            psr = sciCube.hdr['PIXAR_SR']*u.sr
            psrScale = psr / u.pixel
            psa = np.sqrt(psr.to(u.arcsec**2))
            
            #Build the spectral axis
            wv0 = sciCube.hdr['CRVAL3']
            dwv = sciCube.hdr['CDELT3']

            #Generate the wavelength array
            wvls = np.arange(dnpts)*dwv + wv0

            print('Center coordinates are ({:d},{:d})'.format(sciCube.xcenter,sciCube.ycenter))

            #Convert annulus width to pixels
            pixscale = u.pixel_scale(psa/u.pixel)
            innerRadPix = (innerRadius).to(u.pixel,pixscale)
            outerRadPix = (outerRadius).to(u.pixel,pixscale)

            #Define the aperture
            apCen = (sciCube.xcenter,sciCube.ycenter)
            if innerRadPix == 0:
                apEx = CircularAperture(apCen, r = outerRadPix.value)
            else:
                apEx = CircularAnnulus(apCen, r_in = innerRadPix.value, r_out = outerRadPix.value)

            #Perform the extraction, one spectral pixel at a time
            spec = np.zeros(dnpts)
            sigma = np.zeros(dnpts)

            for i in range(dnpts):
                sci_img = data[i,:,:]*u.MJy/u.sr
                err_img = derr[i,:,:]*u.MJy/u.sr
                apPhot = ApertureStats(data = sci_img, aperture = apEx, error = err_img, sum_method='subpixel', subpixels=5)

                flux_jy = (apPhot.mean*psrScale).to(u.Jy/u.pixel)
                noise_jy = (apPhot.std*psrScale/np.sqrt(apPhot.sum_aper_area.value)).to(u.Jy/u.pixel)

                spec[i] = flux_jy.value
                sigma[i] = noise_jy.value

            if withPlots:
                fig, axes = plt.subplots(1,2,figsize=(20,10))
                fig.subplots_adjust(hspace=0.45,wspace=0.15)
                axes[0].imshow(sciCube.cdata,origin='lower',cmap='viridis',interpolation='none')
                axes[0].plot(sciCube.xcenter,sciCube.ycenter,marker='+',markersize=8,color='r')
                axes[0].set_title('Extraction Region for Cube #{}'.format(j))
                apEx.plot(ax=axes[0], color='red',lw=2)

                axes[1].plot(wvls,spec)
                axes[1].set_xlabel('Wavelength ($\mu$m)')
                axes[1].set_ylabel('Flux (Jy)')
                axes[1].set_title('Extracted Spectrum for Cube #{}'.format(j))
                axes[1].set_xlim(waveLo.value,waveUp.value)
                plt.show()
            print(len(wvls),len(spec))

            combined_waves.append(wvls.tolist())
            combined_specs.append(spec.tolist())
            combined_sigmas.append(sigma.tolist())

        if len(cubeFile) == 1:
            wvls, spec, sigma = combined_waves[0], combined_specs[0], combined_sigmas[0]
        else:
            wvls, spec, sigma = subchannel_splice(combined_waves, combined_specs, combined_sigmas)

        #Only extract the wavelength region of interest
        wv_region = np.where((wvls>waveLo.value) & (wvls<waveUp.value))

        #Save the file
        specFile = 'JWST-Extract-{:.2f}-arcsecInnerRadius-{:.2f}-arcsecOuterRadius-{:.2f}um-to-{:.2f}um.txt'.format(innerRadius.value,outerRadius.value,waveLo.value,waveUp.value)
        #Get header observation information
        obsInfo = readHeader(cubeFile[0])
        with open(specFile, 'w') as fn:
            #Create headers with extract information
            fn.write('#Target name {}\n'.format(obsInfo.target))
            fn.write('#Obs. Start {} {}\n'.format(obsInfo.dateBeg,obsInfo.timeBeg))
            fn.write('#Obs. End {} {}\n'.format(obsInfo.dateEnd,obsInfo.timeEnd))
            fn.write('#Lower wavelength (um) {}\n'.format(waveLo.value))
            fn.write('#Upper wavelength (um) {}\n'.format(waveUp.value))
            fn.write('#Center pixel for extract (x,y) = ({},{})\n'.format(sciCube.xcenter,sciCube.ycenter))
            fn.write('#Pixel scale (arcsec/pixel) {}\n'.format(psa.value))
            fn.write('#Inner annulus radius (arcsec) {}\n'.format(innerRadius.value))
            fn.write('#Outer annulus radius (arcsec) {}\n'.format(outerRadius.value))
            fn.write('#Wave (micron) Flux (Jy) Noise (Jy)\n')

            for w, s, e in zip(wvls[wv_region],spec[wv_region],sigma[wv_region]):
                fn.write('{} {} {}\n'.format(w,s,e))

        #Plot the aperture if desired
        if withPlots:
            fig, axes = plt.subplots(1,1,figsize=(10,10))

            axes.plot(wvls[wv_region],spec[wv_region])
            axes.set_xlabel('Wavelength ($\mu$m)')
            axes.set_ylabel('Flux (Jy)')
            axes.set_title('Spliced Extracted Spectrum')
            plt.show()