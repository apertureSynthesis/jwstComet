import os,sys
import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
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
                    specFile = specStem+'-{:.2f}-arcsecRadAp-{:.1f}-arcsecXoff-{:.1f}-arcsecYoff-{:.2f}um-to-{:.2f}um.txt'.format(radAp.value,dxArc.value,dyArc.value,waveLo.value,waveUp.value)
                    resFile = specFile[:-3]+'.retrieval-results.txt'
                    beam = Beam()
                    beamExtract = beam.extractSpec(cubeFiles=cubeFiles, specFile=specFile, waveLos=waveLos, waveUps=waveUps, radAp=radAp, xOffset=dxArc, yOffset=dyArc, smooth=smooth)
                    beamModel = runPSG()
                    beamModel.getModels(specFile=specFile, resFile=resFile, name=name, objectType=objectType, composition=composition, retrieval=retrieval, mode='mapping', withPlots=True, key=key)
                    
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




