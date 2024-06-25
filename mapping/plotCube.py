import os,sys
import numpy as np
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from matplotlib import pyplot as plt
import pandas as pd
from glob import glob
import astropy.units as u
from jwstComet import extraction, utils

class plotCube(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    
    def polynomial_fitting(self,wave,spec,order=2,iters=5,nsigma=2.5,withPlots=False):
        """
        Function to remove continuum with a polynomial fit
        """
        fit = fitting.LevMarLSQFitter()
        m = models.Polynomial1D(degree=order)
        or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=iters, sigma=nsigma)

        f, mask = or_fit(m, wave, spec)
        filtered_data = np.ma.masked_array(spec, mask=~mask)

        if withPlots:
            fig,axes = plt.subplots(2,1)
            axes[0].plot(wave,spec)
            axes[0].plot(wave,filtered_data,"r.")
            axes[0].plot(wave,f(wave))
            axes[1].plot(wave,spec-f(wave))
            plt.show()

        return spec - f(wave)
    
    def makePlots(self,cubeFiles,specStem,csvFile,waveLo,waveUp,radAp,withPlots=False):
        """
        Read in a JWST IFU cube. Find the photocenter. Extract spectra across the
        entire cube. Subtract continuum. Plot the results.
        """

        #Read in the first cube to serve as a coordinate reference
        sciCube = utils.readCube(cubeFiles[0])

        retrieval_x_indexes   = []
        retrieval_x_offsets   = []
        retrieval_y_indexes   = []
        retrieval_y_offsets   = []
        extracted_spectra = []
        subtracted_spectra = []

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
        for x in range(0, sciCube.xs):
            for y in range(0, sciCube.ys):
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
                    beam = extraction.Beam()
                    beamExtract = beam.extractSpec(cubeFiles=cubeFiles, specFile=specFile, waveLo=waveLo, waveUp=waveUp, radAp=radAp, xOffset=dxArc, yOffset=dyArc)

                    #Perform the continuum subtraction
                    wave, spec, err = np.loadtxt(specFile, unpack=1)
                    if np.isnan(np.sum(spec)):
                        continue

                    spec_clipped = sigma_clip(spec,sigma=10,maxiters=5)
                    spec_sub = self.polynomial_fitting(wave, spec_clipped, withPlots=withPlots)

                    #Add the position to the index list
                    retrieval_x_indexes.append(x)
                    retrieval_y_indexes.append(y)
                    retrieval_x_offsets.append(dxArc.value)
                    retrieval_y_offsets.append(dyArc.value)
                    extracted_spectra.append(np.nansum(spec))
                    subtracted_spectra.append(np.nansum(spec_sub))

        #Save the results to a CSV                  
        #Create a dataframe to store
        df = pd.DataFrame()

        df['X-Index'] = retrieval_x_indexes
        df['Y-Index'] = retrieval_y_indexes
        df['X-Offset'] = retrieval_x_offsets
        df['Y-Offset'] = retrieval_y_offsets
        df['Spectra'] = extracted_spectra
        df['Sub-spectra'] = subtracted_spectra

        df.to_csv(csvFile,index=False)      

    def plotMaps(self,csvFile):
        """
        Plot out the results of a map for each retrieved value
        """

        #Read in the saved results
        df = pd.read_csv(csvFile)


        #Split off the retrieved variables and sigmas
        df_retrieved = df.drop(['X-Index','Y-Index','X-Offset','Y-Offset'], axis=1)
        df_values = df_retrieved.iloc[:,:]

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



