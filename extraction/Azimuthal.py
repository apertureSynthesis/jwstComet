import numpy as np
import astropy.units as u
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats, aperture_photometry
from jwstComet.utils import readCube, subchannel_splice, readHeader

class Azimuthal(object):

    def __init__(self):
        super().__init__()
        self.name = self.__class__.__name__

    @u.quantity_input(innerRadius = u.arcsec, outerRadius = u.arcsec)
    def extractSpec(self,cubeFiles,specFile,waveLo,waveUp,innerRadius=0.0*u.arcsec,outerRadius=0.1*u.arcsec,
                    method='subpixel',subpixels=10,mask=None,withPlots=False,split=None):
        """
        Extract a spectrum within a specified annulus, take the azimuthal average, and save it to a text file.
        Plot the aperture and extracted spectrum if desired.
        
        Inputs
            cubeFiles - array of file paths pointing to the *s3d.fits datacube files from which we want to extract a spectrum
            specFile - name of the output file containing the extracted spectrum
            waveLo - lowest wavelength for extraction. preferred unit is microns. can be a list or single value
            waveUp - highest wavelength for extraction. preferred unit is microns. can be a list or single value
            innerRadius - inner radius (arcsec) for the annulus
            outerRadius - outer radius (arcsec) for the annulus
            method - extraction method for photutils, choosing from options of 'subpixel', 'center', or 'exact'
            subpixels - number of subpixels to be used if subpixel extraction is chosen
            mask - dictionary containing lower and upper wavelength range to be masked. optional
            withPlots - whether we plot the results. optional
            split - whether we split the array and only analyze a subsection. choices are 'upper', 'lower', 'left', 'right'. optional

        Outputs
            Saves an ASCII file containing the extracted spectrum and header information. Optionally shows (but does not save) plots.
        """
        if type(waveLo) is not list: 
            waveLo = [waveLo]
        if type(waveUp) is not list: 
            waveUp = [waveUp]

        combined_waves = []
        combined_specs = []
        combined_sigmas = []

        for j in range(len(cubeFiles)):

            #Read in the file and get relevant information and data
            sciCube = readCube(cubeFiles[j])
            data = sciCube.data
            derr = sciCube.derr
            cdata = sciCube.cdata
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

            #Check whether we are splitting the chip
            if split == 'left':
                data[:,:,sciCube.xcenter:] = np.nan
                derr[:,:,sciCube.xcenter:] = np.nan
                cdata[:,sciCube.xcenter:] = np.nan
            elif split == 'right':
                data[:,:,:sciCube.xcenter] = np.nan
                derr[:,:,:sciCube.xcenter] = np.nan
                cdata[:,:sciCube.xcenter] = np.nan
            elif split == 'lower':
                data[:,sciCube.ycenter:,:] = np.nan
                derr[:,sciCube.ycenter:,:] = np.nan
                cdata[sciCube.ycenter:,:] = np.nan
            elif split == 'upper':
                data[:,:sciCube.ycenter,:] = np.nan
                derr[:,:sciCube.ycenter,:] = np.nan
                cdata[:sciCube.ycenter,:] = np.nan

            for i in range(dnpts):
                sci_img = data[i,:,:]*u.MJy/u.sr * psr
                err_img = derr[i,:,:]*u.MJy/u.sr * psr
                wgt_img = 1./(err_img**2)
                #Take the weighted average: Sum(w_i * x_i) / Sum(w_i), w_i = 1/noise_i^2
                #Uncertainty in weighted average: 1 / sqrt(Sum(w_i))
                #Taylor, An Introduction to Error Analysis
                if method == 'subpixel':
                    apPhot = ApertureStats(data = sci_img*wgt_img, aperture = apEx, sum_method='subpixel', subpixels=subpixels)
                    wtPhot = ApertureStats(data = wgt_img, aperture = apEx, sum_method='subpixel', subpixels=subpixels)
                elif method == 'center':
                    apPhot = ApertureStats(data = sci_img*wgt_img, aperture = apEx, sum_method='center')
                    wtPhot = ApertureStats(data = wgt_img, aperture = apEx, sum_method='center')
                elif method == 'exact':
                    apPhot = ApertureStats(data = sci_img*wgt_img, aperture = apEx, sum_method='exact')
                    wtPhot = ApertureStats(data = wgt_img, aperture = apEx, sum_method='exact')                    


                flux_jy = ((apPhot.sum / wtPhot.sum)).to(u.Jy)
                noise_jy = ((1./np.sqrt(wtPhot.sum))).to(u.Jy)

                #If the annulus is larger than one pixel, scale the fluxes by the width
                #to compensate for the radial decay in flux with distance
                if (outerRadPix - innerRadPix).value > 1:
                    annulusWidth = (outerRadPix - innerRadPix).value
                    flux_jy *= annulusWidth
                    noise_jy *= annulusWidth

                spec[i] = flux_jy.value
                sigma[i] = noise_jy.value

            if withPlots:
                fig, axes = plt.subplots(1,2)
                fig.subplots_adjust(hspace=0.45,wspace=0.45)
                axes[0].imshow(cdata,origin='lower',cmap='viridis',interpolation='none')
                axes[0].plot(sciCube.xcenter,sciCube.ycenter,marker='+',markersize=6,color='r')
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
        wv_region = []
        for waveLo, waveUp in zip(waveLo,waveUp):
            wvr = np.where((wvls>waveLo.value) & (wvls<waveUp.value))
            wv_region = np.concatenate((wv_region,wvr[0]))

        wv_region = wv_region.astype(int)

        #Mask pixels with very large sigma if requested
        if mask != None:
            mask_keys = list(mask.keys())
            for key in range(len(mask_keys)):
                maskLo = mask[mask_keys[key]]['start']
                maskUp = mask[mask_keys[key]]['stop']
                mask_region = np.where((wvls>maskLo.value) & (wvls<maskUp.value))
                sigma[mask_region] = np.nanmax(spec[mask_region])


        #Save the file
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
            fn.write('#Center pixel for extract (x,y) = ({},{})\n'.format(sciCube.xcenter,sciCube.ycenter))
            fn.write('#Pixel scale (arcsec/pixel) {}\n'.format(psa.value))
            fn.write('#Inner annulus radius (arcsec) {}\n'.format(innerRadius.value))
            fn.write('#Outer annulus radius (arcsec) {}\n'.format(outerRadius.value))
            fn.write('#Wave (micron) Flux (Jy) Noise (Jy)\n')

            for w, s, e in zip(wvls[wv_region],spec[wv_region],sigma[wv_region]):
                fn.write('{} {} {}\n'.format(w,s,e))

        #Plot the aperture if desired
        if withPlots:
            fig, axes = plt.subplots(1,1)

            axes.errorbar(wvls[wv_region],spec[wv_region],sigma[wv_region],ecolor='r')
            axes.set_xlabel('Wavelength ($\mu$m)')
            axes.set_ylabel('Flux (Jy)')
            axes.set_title('Spliced Extracted Spectrum')
            plt.show()