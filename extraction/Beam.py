
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.convolution import convolve, Box2DKernel
from photutils.aperture import CircularAperture, aperture_photometry
from jwstComet.utils import readCube, subchannel_splice, readHeader

class Beam(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    @u.quantity_input(radAp=u.arcsec, xOffset=u.arcsec, yOffset=u.arcsec)
    def extractSpec(self,cubeFiles,specFile,waveLos,waveUps,radAp,xOffset=0.0*u.arcsec,yOffset=0.0*u.arcsec,smooth=None,mask=None,withPlots=False):
        """
        Extract a spectrum at the specified position and save it to a text file.
        Plot the aperture and extracted spectrum if desired
        """
        if type(waveLos) is not list:
            waveLos = [waveLos]
        
        if type(waveUps) is not list:
            waveUps = [waveUps]

        # try:
        #     if waveLo.value >= waveUp.value:
        #         raise ValueError('Lower wavelength must be smaller than upper wavelength')
        # except (ValueError):
        #     exit('Update wavelength range')

        combined_waves = []
        combined_specs = []
        combined_sigmas = []

        for j in range(len(cubeFiles)):
            #Read in the file and get relevant information and data
            sciCube = readCube(cubeFiles[j])
            #Data
            data = sciCube.data
            #Noise
            derr = sciCube.derr
            #Collapsed cube (for plotting)
            cdata = sciCube.cdata
            #Length of spectral axis
            dnpts = data.shape[0]

            if smooth != None:
                cdata = convolve(sciCube.cdata,Box2DKernel(smooth))

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
                if smooth != None:
                    sci_img = convolve(data[i,:,:],Box2DKernel(smooth))*u.MJy/u.sr
                    err_img = convolve(data[i,:,:],Box2DKernel(smooth))*u.MJy/u.sr
                else:
                    sci_img = data[i,:,:]*u.MJy/u.sr
                    err_img = derr[i,:,:]*u.MJy/u.sr
                apPhot = aperture_photometry(data = sci_img, apertures = apEx, error = err_img, method='subpixel', subpixels=5)

                flux_jy = (apPhot['aperture_sum'][0]*psrScale).to(u.Jy/u.pixel)
                noise_jy = (apPhot['aperture_sum_err'][0]*psrScale).to(u.Jy/u.pixel)

                spec[i] = flux_jy.value
                sigma[i] = noise_jy.value

            if withPlots:
                fig, axes = plt.subplots(1,2)
                fig.subplots_adjust(hspace=0.45,wspace=0.45)
                axes[0].imshow(cdata,origin='lower',cmap='viridis',interpolation='none')
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

        if len(cubeFiles) == 1:
            wvls, spec, sigma = np.array(combined_waves[0]), np.array(combined_specs[0]), np.array(combined_sigmas[0])
        else:
            wvls, spec, sigma = subchannel_splice(combined_waves, combined_specs, combined_sigmas)



        #Only extract the wavelength region of interest
        #wv_region = np.where((wvls>waveLo.value) & (wvls<waveUp.value))
        wv_region = []
        for waveLo, waveUp in zip(waveLos,waveUps):
            wvr = np.where((wvls>waveLo.value) & (wvls<waveUp.value))
            wv_region = np.concatenate((wv_region,wvr))

        wv_region = wv_region.astype(int)

        #Mask pixels with very large sigma if requested
        if mask != None:
            mask_keys = list(mask.keys())
            for key in range(len(mask_keys)):
                maskLo = mask[mask_keys[key]]['start']
                maskUp = mask[mask_keys[key]]['stop']
                mask_region = np.where((wvls>maskLo.value) & (wvls<maskUp.value))
                sigma[mask_region] *= 1e25

        #Save to a file
        #Get header observation information
        obsInfo = readHeader(cubeFiles[0])
        with open(specFile, 'w') as fn:
            #Create headers with extract information
            fn.write('#Target name {}\n'.format(obsInfo.target))
            fn.write('#Obs. Start {} {}\n'.format(obsInfo.dateBeg,obsInfo.timeBeg))
            fn.write('#Obs. End {} {}\n'.format(obsInfo.dateEnd,obsInfo.timeEnd))
            fn.write('#Instrument used {}\n'.format(obsInfo.instrument))
            fn.write('#Grating used {}\n'.format(obsInfo.setting))
            fn.write('#Lower wavelength (um) {}\n'.format(waveLo.value))
            fn.write('#Upper wavelength (um) {}\n'.format(waveUp.value))
            fn.write('#Spectral plate scale (um) {}\n'.format(dwv))
            fn.write('#Center pixel for extract (x,y) = {},{}\n'.format(sciCube.xcenter,sciCube.ycenter))
            fn.write('#Pixel scale (arcsec/pixel) {}\n'.format(psa.value))
            fn.write('#Aperture radius (arcsec) {}\n'.format(radAp.value))
            fn.write('#X offset (arcsec) {}\n'.format(xOffset.value))
            fn.write('#Y offset (arcsec) {}\n'.format(yOffset.value))
            fn.write('#Wave (micron) Flux (Jy) Noise (Jy)\n')

            for w, s, e in zip(wvls[wv_region],spec[wv_region],sigma[wv_region]):
                fn.write('{} {} {}\n'.format(w,s,e))

        #Plot the aperture if desired
        if withPlots:
            fig, axes = plt.subplots(1,1)

            axes.errorbar(wvls[wv_region],spec[wv_region],sigma[wv_region],ecolor='r')
            axes.set_xlabel('Wavelength ($\mu$m)')
            axes.set_ylabel('Flux (Jy)')
            axes.set_title('Extracted Spectrum for {:.2f} radius aperture\n xOffset = {:.2f}, yOffset = {:.2f}'.format(radAp,xOffset,yOffset))
            plt.show()