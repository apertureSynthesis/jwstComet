import os,sys
import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box2DKernel
from photutils.aperture import RectangularAperture, aperture_photometry
from jwstComet.extraction import Beam
from jwstComet.modeling import runPSG, readPSG
from jwstComet.utils import readCube

class Mapping(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @u.quantity_input(radAp=u.arcsec)
    def makeMaps(self,cubeFiles,specStem,csvFile,waveLos,waveUps,radAp,name,objectType,composition,retrieval,smooth=None,box=None,key=None):
        """
        Read in a JWST IFU cube. Find the photocenter. Extract spectra across the
        entire cube. Send them to the PSG for analysis. Plot the model results and extracted spectrum.
        Save the results to a CSV file.
        """

        #Read in the first cube to serve as a coordinate reference
        sciCube = readCube(cubeFiles[0])

        #Crop a portion if desired
        if box != None:
            x0 = sciCube.xcenter - box
            xf = sciCube.xcenter + box
            y0 = sciCube.ycenter - box
            yf = sciCube.ycenter + box
        else:
            x0 = 0
            xf = sciCube.xs
            y0 = 0
            yf = sciCube.ys

        retrieval_x_indexes   = []
        retrieval_x_offsets   = []
        retrieval_y_indexes   = []
        retrieval_y_offsets   = []
        retrieval_variables   = []
        retrieval_values      = []
        retrieval_sigmas      = []

        #Work out how much we are spatially binning the results
        #Solid angle subtended by a square pixel (steradians)
        psr = sciCube.hdr['PIXAR_SR']*u.sr
        #Convert to pixel side length in arcseconds
        psa = np.sqrt(psr.to(u.arcsec**2))
        #Assign a pixel scale for conversion between pixels and arcseconds
        pixScale = u.pixel_scale(psa/u.pixel)
        #Bin factor
        bin_factor =  int(np.round(radAp / psa))


        #Now extract the spatially binned (if requested) maps
        for x in range(x0, xf):
            for y in range(y0, yf):
                #Check whether we are on the chip
                if sciCube.wmap[50,y,x] != 0:
                    #Calculate x_offset and y_offset in arcseconds
                    dxPix = (x - sciCube.xcenter)*u.pixel
                    dxArc = dxPix.to(u.arcsec,pixScale)
                    dyPix = (y - sciCube.ycenter)*u.pixel
                    dyArc = dyPix.to(u.arcsec,pixScale)


                    #Perform the extract
                    specFile = specStem+'-{:.2f}-arcsecRadAp-{:.1f}-arcsecXoff-{:.1f}-arcsecYoff-{:.2f}um-to-{:.2f}um.txt'.format(radAp.value,dxArc.value,dyArc.value,min(waveLos).value,max(waveUps).value)
                    resFile = specFile[:-3]+'.retrieval-results.txt'
                    beam = Beam()
                    beamExtract = beam.extractSpec(cubeFiles=cubeFiles, specFile=specFile, waveLos=waveLos, waveUps=waveUps, radAp=radAp, xOffset=dxArc, yOffset=dyArc, smooth=smooth)
                    beamModel = runPSG()
                    beamModel.getModels(specFile=specFile, resFile=resFile, name=name, objectType=objectType, composition=composition, retrieval=retrieval, mode='beam', withPlots=True, key=key)
                    
                    try:
                        results = readPSG(resFile)
                        retrieval_variables.append(results.retrieval_variables)
                        retrieval_values.append(results.retrieval_values)
                        retrieval_sigmas.append(results.retrieval_sigmas)

                        #Add the position to the index list
                        retrieval_x_indexes.append(x)
                        retrieval_y_indexes.append(y)
                        retrieval_x_offsets.append(dxArc.value)
                        retrieval_y_offsets.append(dyArc.value)

                        #Save results to a file, or update the existing file
                        if os.path.exists(csvFile):
                            df1 = pd.read_csv(csvFile)

                            df = pd.DataFrame()
                            df['X-Index'] = [x]
                            df['Y-Index'] = [y]
                            df['X-Offset'] = [dxArc.value]
                            df['Y-Offset'] = [dyArc.value]

                            for i in range(len(results.retrieval_variables)):
                                df[results.retrieval_variables[i]] = [results.retrieval_values[i]]
                                df['sigma-'+results.retrieval_variables[i]] = [results.retrieval_sigmas[i]]

                            df2 = pd.concat([df1,df])
                            os.system('rm {}'.format(csvFile))
                            df2.to_csv(csvFile,index=False)

                        else:
                            df = pd.DataFrame()
                            df['X-Index'] = [x]
                            df['Y-Index'] = [y]
                            df['X-Offset'] = [dxArc.value]
                            df['Y-Offset'] = [dyArc.value]

                            for i in range(len(results.retrieval_variables)):
                                df[results.retrieval_variables[i]] = [results.retrieval_values[i]]
                                df['sigma-'+results.retrieval_variables[i]] = [results.retrieval_sigmas[i]]

                            df.to_csv(csvFile,index=False)

                    except:
                        pass

    def makeImage(self,cubeFiles,specStem,csvFile,waveLos,waveUps,smooth=None,box=None,waveContLos=None,waveContUps=None):
        """
        Plot an integrated intensity map. Subtract continuum if desired.
        Save results to a CSV and PDF
        """

        if type(waveLos) is not list:
            waveLos = [waveLos]

        if type(waveUps) is not list:
            waveUps = [waveUps]

        if (waveContLos is not None) & (type(waveContLos) is not list):
            waveContLos = [waveContLos]

        if (waveContLos is not None) & (type(waveContLos) is not list):
            waveContLos = [waveContLos]

        #Read in the first cube to serve as a coordinate reference
        sciCube = readCube(cubeFiles[0])

        #Crop a portion if desired
        if box != None:
            x0 = sciCube.xcenter - box
            xf = sciCube.xcenter + box
            y0 = sciCube.ycenter - box
            yf = sciCube.ycenter + box
        else:
            x0 = 0
            xf = sciCube.xs
            y0 = 0
            yf = sciCube.ys

        retrieval_x_indexes   = []
        retrieval_x_offsets   = []
        retrieval_y_indexes   = []
        retrieval_y_offsets   = []
        retrieval_values      = []
        retrieval_sigmas      = []

        #Work out how much we are spatially binning the results
        #Solid angle subtended by a square pixel (steradians)
        psr = sciCube.hdr['PIXAR_SR']*u.sr
        #Convert to pixel side length in arcseconds
        psa = np.sqrt(psr.to(u.arcsec**2))
        #Assign a pixel scale for conversion between pixels and arcseconds
        pixScale = u.pixel_scale(psa/u.pixel)

        #Build the spectral axis
        wv0 = sciCube.hdr['CRVAL3']
        dwv = sciCube.hdr['CDELT3']

        #Generate the wavelength array
        dnpts = sciCube.data.shape[0]
        wvls = np.arange(dnpts)*dwv + wv0

        #smooth if desired
        if smooth != None:
            cdata = convolve(sciCube.cdata,Box2DKernel(smooth))
           
        #Now extract the spatially binned (if requested) maps
        for x in range(x0, xf):
            for y in range(y0, yf):
                #Check whether we are on the chip
                if sciCube.wmap[50,y,x] != 0:
                    #Calculate x_offset and y_offset in arcseconds
                    dxPix = (x - sciCube.xcenter)*u.pixel
                    dxArc = dxPix.to(u.arcsec,pixScale)
                    dyPix = (y - sciCube.ycenter)*u.pixel
                    dyArc = dyPix.to(u.arcsec,pixScale)

                    apEx = RectangularAperture((y,x),w=1,h=1)

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


                    #Subtract continuum
                    if (waveContLos is not None) & (waveContUps is not None):
                        cont_region = []
                        for waveContLo, waveContUp in zip(waveContLos, waveContUps):
                            cwv = np.where((wvls>waveContLo.value) & (wvls<waveContUp.value))
                            cont_region = np.concatenate((cont_region,cwv[0]))
                        
                        #Find the median value of the continuum
                        cont_med = np.nanmedian(spec[cont_region])

                        spec -= cont_med

                    #Save to a file
                    #Perform the extract
                    specFile = specStem+'-{:.1f}-arcsecXoff-{:.1f}-arcsecYoff-{:.2f}um-to-{:.2f}um.txt'.format(dxArc.value,dyArc.value,min(waveLos).value,max(waveUps).value)

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



    def plotMaps(self,csvFile):
        """
        Plot out the results of a map for each retrieved value
        """

        #Read in the saved results
        df = pd.read_csv(csvFile)


        #Split off the retrieved variables and sigmas
        df_retrieved = df.drop(['X-Index','Y-Index','X-Offset','Y-Offset'], axis=1)
        df_values = df_retrieved.iloc[:,::2]

        #Set up the plots
        fig, axes = plt.subplots(1,len(df_values.keys()),figsize=(10,10))
        fig.subplots_adjust(hspace=0.25,wspace=0.01)

        if (len(df_values.keys()) > 1):
            axes = axes.ravel()

            #Now plot out one at a time:
            for i in range(len(df_values.keys())):
                plot_array = np.zeros(((df['Y-Index'].max()+1),(df['X-Index'].max()+1)))
                for x,y,value in zip(df['X-Index'],df['Y-Index'],df_values[df_values.columns[i]]):
                    plot_array[y,x] = value
                im = axes[i].imshow(plot_array,origin='lower',cmap='viridis')
                axes[i].set_title(df_values.columns[i])
                plt.colorbar(im, ax=axes[i])

        else:
                plot_array = np.zeros(((df['X-Index'].max()+1),(df['Y-Index'].max()+1)))
                for x,y,value in zip(df['X-Index'],df['Y-Index'],df_values[df_values.columns[0]]):
                    plot_array[y,x] = value
                im = axes.imshow(plot_array,origin='lower',cmap='viridis')
                axes.set_title(df_values.columns[0])
                plt.colorbar(im, ax = axes)




