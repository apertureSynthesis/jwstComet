import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAperture, aperture_photometry
from jwstComet.utils import readCube, subchannel_splice, readHeader

class beam(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    @u.quantity_input(waveLo=u.um, waveUp=u.um, radAp=u.arcsec, xOffset=u.arcsec, yOffset=u.arcsec)
    def extractSpec(self,cubeFile,waveLo,waveUp,radAp,xOffset=0.0*u.arcsec,yOffset=0.0*u.arcsec,withPlots=False):
        """
        Extract a spectrum at the specified position and save it to a text file.
        Plot the aperture and extracted spectrum if desired
        """

        try:
            if waveLo.value >= waveUp.value:
                raise ValueError('Lower wavelength must be smaller than upper wavelength')
        except (ValueError):
            exit('Update wavelength range')

        combined_waves = []
        combined_specs = []
        combined_sigmas = []

        for j in range(len(cubeFile)):
            #Read in the file and get relevant information and data
            sciCube = readCube(cubeFile[j])
            #Data
            data = sciCube.data
            #Noise
            derr = sciCube.derr
            #Length of spectral axis
            dnpts = data.shape[0]

            #Pixel scale - convert from steradians to arcsec**2
            psr = sciCube.hdr['PIXAR_SR']*u.sr
            psrScale = psr / u.pixel
            psa = np.sqrt(psr.to(u.arcsec**2))
            
            #Build the spectral axis
            wv0 = sciCube.hdr['CRVAL3']
            dwv = sciCube.hdr['CDELT3']

            #Generate the wavelength array
            wvls = np.arange(dnpts)*dwv + wv0

            print('Center coordinates are ({:d},{:d})'.format(sciCube.xcenter,sciCube.ycenter))

            #Convert aperture radius, horizontal and vertical offsets from arcseconds to pixels
            pixscale = u.pixel_scale(psa/u.pixel)
            radpix = (radAp).to(u.pixel,pixscale)
            xOffpix  = (xOffset).to(u.pixel,pixscale)
            yOffpix  = (yOffset).to(u.pixel,pixscale)

            #Define the aperture
            apCen = (sciCube.xcenter+xOffpix.value, sciCube.ycenter+yOffpix.value)
            apEx  = CircularAperture(apCen, r = radpix.value)

            #Perform the extraction, one spectral pixel at a time
            spec = np.zeros(dnpts)
            sigma = np.zeros(dnpts)

            for i in range(dnpts):
                sci_img = data[i,:,:]*u.MJy/u.sr
                err_img = derr[i,:,:]*u.MJy/u.sr
                apPhot = aperture_photometry(data = sci_img, apertures = apEx, error = err_img, method='subpixel', subpixels=5)

                flux_jy = (apPhot['aperture_sum'][0]*psrScale).to(u.Jy/u.pixel)
                noise_jy = (apPhot['aperture_sum_err'][0]*psrScale).to(u.Jy/u.pixel)

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
                plt.show()

            combined_waves.append(wvls.tolist())
            combined_specs.append(spec.tolist())
            combined_sigmas.append(sigma.tolist())

        if len(cubeFile) == 1:
            wvls, spec, sigma = np.array(combined_waves[0]), np.array(combined_specs[0]), np.array(combined_sigmas[0])
        else:
            wvls, spec, sigma = subchannel_splice(combined_waves, combined_specs, combined_sigmas)



        #Only extract the wavelength region of interest
        wv_region = np.where((wvls>waveLo.value) & (wvls<waveUp.value))

        #Save to a file
        specFile = 'JWST-Extract-{:.2f}-arcsecRadAp-{:.1f}-arcsecXoff-{:.1f}-arcsecYoff-{:.2f}um-to-{:.2f}um.txt'.format(radAp.value,xOffset.value,yOffset.value,waveLo.value,waveUp.value)
        #Get header observation information
        obsInfo = readHeader(cubeFile[0])
        with open(specFile, 'w') as fn:
            #Create headers with extract information
            fn.write('#Target name {}\n'.format(obsInfo.target))
            fn.write('#Obs. Start {} {}\n'.format(obsInfo.dateBeg,obsInfo.timeBeg))
            fn.write('#Obs. End {} {}\n'.format(obsInfo.dateEnd,obsInfo.timeEnd))
            fn.write('#Lower wavelength (um) {}\n'.format(waveLo.value))
            fn.write('#Upper wavelength (um) {}\n'.format(waveUp.value))
            fn.write('Center pixel for extract (x,y) = {},{}\n'.format(sciCube.xcenter,sciCube.ycenter))
            fn.write('#Aperture radius (arcsec) {}\n'.format(radAp.value))
            fn.write('#X offset (arcsec) {}\n'.format(xOffset.value))
            fn.write('#Y offset (arcsec) {}\n'.format(yOffset.value))
            fn.write('#Wave (micron) Flux (Jy) Noise (Jy)\n')

            for w, s, e in zip(wvls[wv_region],spec[wv_region],sigma[wv_region]):
                fn.write('{} {} {}\n'.format(w,s,e))

        #Plot the aperture if desired
        if (withPlots and len(cubeFile) > 1):
            fig, axes = plt.subplots(1,1,figsize=(10,10))

            axes.plot(wvls[wv_region],spec[wv_region])
            axes.set_xlabel('Wavelength ($\mu$m)')
            axes.set_ylabel('Flux (Jy)')
            axes.set_title('Spliced Extracted Spectrum')
            plt.show()