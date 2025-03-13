import os,sys
import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from astropy.convolution import convolve, Box2DKernel
from photutils.aperture import RectangularAperture, aperture_photometry
from jwstComet.extraction import Beam
from jwstComet.modeling import runPSG, readPSG
from jwstComet.utils import readCube

class Mapping(object):

    def __init__(self):
        super().__init__()
        self.name = self.__class__.__name__

    @u.quantity_input(radAp=u.arcsec)
    def makeMaps(self,cubeFiles,specStem,csvFile,waveLo,waveUp,radAp,name,objectType,composition,retrieval,withCont=False,smooth=None,box=None,withEph=True,local=True,tempFix=False,tempFile=None):
        """
        Read in a JWST IFU cube. Find the photocenter. Extract spectra across the
        entire cube. Send them to the PSG for analysis. Plot the model and extracted spectrum.
        Save the results to a CSV file.

        Inputs
            cubeFiles - array of file paths pointing to the *s3d.fits datacube files from which we want to extract a spectrum
            specStem - stem for extracted spectra output file names
            csvFile - name of CSV file for saving pixel-by-pixel results
            waveLo - lowest wavelength for extraction. preferred unit is microns. can be a list or single value
            waveUp - highest wavelength for extraction. preferred unit is microns. can be a list or single value
            radAp - radius of the circular extraction aperture (arcsec) or [x,y] lengths of the rectangular extraction aperture (arcsec,arcsec)
            name - name of the comet or asteroid
            objectType - type of small body: comet or asteroid
            composition - dictionary containing compositional information for building the PSG model atmosphere
            retrieval - dictionary containing quantities to be retrieved for each PSG model run
            withCont - whether we are asking the PSG to simulate the continuum or instead simply subtract a baseline
            smooth - kernel length for Box2DKernel smoothing of the cube. optional.
            box - length of a box to crop the image around the photocenter. optional
            withEph - whether we are asking the PSG to retrieve ephemeris parameters or are instead using a local copy.
            local - are we interrogating a local copy of the PSG or instead sending requests to the online server
            tempFix - are we fixing the temperature to pre-determined values? optional
            tempFile - FITS file containing temperature values (as function of pixel index) if using prefixed values

        Outputs
            Saves an ASCII file containing the extracted spectrum and header information at each pixel. Saves the PSG model files from each pixel. 
            Saves a CSV file containing the retrieved values and uncertainties at each pixel. Optionally shows (but does not save) plots.
        """
        if type(waveLo) is not list:
            waveLo = [waveLo]
        
        if type(waveUp) is not list:
            waveUp = [waveUp]

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

        #Solid angle subtended by a square pixel (steradians)
        psr = sciCube.hdr['PIXAR_SR']*u.sr
        #Convert to pixel side length in arcseconds
        psa = np.sqrt(psr.to(u.arcsec**2))
        #Assign a pixel scale for conversion between pixels and arcseconds
        pixScale = u.pixel_scale(psa/u.pixel)

        #If there is a temperature profile provided, read it in and interpolate a model
        if tempFix:
            #Read in a temperature profile stored in a FITS file. Create an interpolated model. May need to change pixel scale - this assumes the temperature profile comes from MIRI/CH1
            dfits = fits.open(tempFile)
            #temperature
            para_temps = dfits[3].data
            #MIRI pixel
            rho = dfits[1].data
            #MIRI pixel scale
            pScale = 0.13*u.arcsec
            #Interpolate
            temp_model = interp1d(rho[:len(para_temps)]*pScale.value,para_temps,kind='linear',fill_value=(para_temps[0],para_temps[-1]),bounds_error=False)

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
                    drArc = np.sqrt(dxArc.value**2 + dyArc.value**2)

                    if tempFix:
                        qtemp = temp_model(drArc)
                        composition['TEMPERATURE']['value'] = qtemp

                    #Perform the extract
                    specFile = specStem+'-{:.2f}-arcsecRadAp-{:.1f}-arcsecXoff-{:.1f}-arcsecYoff-{:.2f}um-to-{:.2f}um.txt'.format(radAp.value,dxArc.value,dyArc.value,min(waveLo).value,max(waveUp).value)
                    resFile = specFile[:-3]+'.retrieval-results.txt'
                    beam = Beam()
                    beamExtract = beam.extractSpec(cubeFiles=cubeFiles, specFile=specFile, waveLo=waveLo, waveUp=waveUp, radAp=radAp, xOffset=dxArc, yOffset=dyArc, mode='rectangle', smooth=smooth, withPlots=True)
                    beamModel = runPSG()
                    beamModel.getModels(specFile=specFile, resFile=resFile, name=name, objectType=objectType, composition=composition, retrieval=retrieval, mode='mapping', withCont=withCont, withPlots=True, withEph=withEph, local=local)
                    
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


    def plotMaps(self,csvFile):
        """
        Plot out the results of a map for each retrieved value.

        Inputs
            csvFile - file containing the retrieved values and sigmas at each pixel to be mapped
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




