import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAperture, aperture_photometry
from jwstComet.utils.readCube import readCube

class beam:

    def __init__(self, cubeFile):
        self.cubeFile = cubeFile
        print(cubeFile)


    @u.quantity_input(waveLo=u.um, waveUp=u.um, radAp=u.arcsec, xOffset=u.arcsec, yOffset=u.arcsec)
    def extractSpec(self,specFile,waveLo,waveUp,radAp,xOffset=0.0*u.arcsec,yOffset=0.0*u.arcsec,withPlots=False):
        """
        Extract a spectrum at the specified position and save it to a text file.
        Plot the aperture and extracted spectrum if desired
        """

        #Read in the file and get relevant information and data
        sciCube = readCube(self.cubeFile)
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

        #Only extract the wavelength region of interest
        wv_region = np.where((wvls>waveLo.value) & (wvls<waveUp.value))

        #Save to a file
        
        with open(specFile, 'w') as fn:
            #Create headers with extract information
            fn.write('#Lower wavelength (um): {}\n'.format(waveLo.value))
            fn.write('#Upper wavelength (um): {}\n'.format(waveUp.value))
            fn.write('#Aperture radius (arcsec): {}\n'.format(radAp.value))
            fn.write('#X offset (arcsec): {}\n'.format(xOffset.value))
            fn.write('#Y offset (arcsec): {}\n'.format(yOffset.value))
            fn.write('Wave (micron) Flux (Jy) Noise (Jy)\n')

            for w, s, e in zip(wvls[wv_region],spec[wv_region],sigma[wv_region]):
                fn.write('{} {} {}\n'.format(w,s,e))

        #Plot the aperture if desired
        if withPlots:
            fig, axes = plt.subplots(1,2,figsize=(10,10))
            axes[0].imshow(sciCube.cdata,origin='lower',cmap='viridis',interpolation='none')
            axes[0].plot(sciCube.xcenter,sciCube.ycenter,marker='+',markersize=8,color='r')
            axes[0].set_title('Extraction Region')
            apEx.plot(ax=axes[0], color='red',lw=2)

            axes[1].plot(wvls[wv_region],spec[wv_region])
            axes[1].set_title('Extracted Spectrum')
            plt.show()