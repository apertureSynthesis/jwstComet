import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
from jwstComet.utils.weightedAverage import weightedAverage

class savePSG(object):
    """
    Given retrieved variables, values, and sigmas, along with geocentric distance and extract location,
    plot out the derived Q-curves for all variables. Save results in a plot and CSV file.
    If desired, convert trace species abundances to Q's and save those to a plot, too.
    """

    
    def saveResults(self, csvFile, retrieval_variables, retrieval_values, retrieval_sigmas, extracts):
        """
        Save PSG results to a CSV file
        """

        #Create dataframe to store
        df = pd.DataFrame()

        #Calculate extract offset positions in km at the comet
        #dist = [(delta).to(u.km)*np.tan(i.to(u.rad)) for i in extracts]
        #extract_dist = [i.value for i in dist]

        df['Distance (arcsec)'] = extracts

        #Now save each value
        for i in range(len(retrieval_variables[0])):
            item = []
            sigma = []
            for j in range(len(retrieval_variables)):
                item.append(retrieval_values[j][i])
                sigma.append(retrieval_sigmas[j][i])


            df[retrieval_variables[0][i]] = item
            df['sigma-'+retrieval_variables[0][i]] = sigma


        df.to_csv(csvFile,index=False)        

    def plotProfiles(self, csvFile, plotFile, plotTitle):

        """
        Plot the native profiles retrieved by PSG for each variable
        """

        #Read in the csv file with the results
        df = pd.read_csv(csvFile)
        
        #Split off the retrieved variables and sigmas
        df_retrieved = df.drop('Distance (arcsec)',axis=1)
        df_values = df_retrieved.iloc[:,::2]
        df_sigmas = df_retrieved.iloc[:,1::2]

        #Set up the plots
        fig, axes = plt.subplots(len(df_values.keys()),1,figsize=(10,10))
        fig.subplots_adjust(hspace=0.25,wspace=0.01)

        if (len(df_values.keys()) > 1):
            axes = axes.ravel()

            #Now plot out one at a time
            for i in range(len(df_values.keys())):
                axes[i].errorbar(df['Distance (arcsec)'][:-1],df_values[df_values.columns[i]][:-1],yerr=df_sigmas[df_sigmas.columns[i]][:-1],color='C0',marker='x',capsize=2,linestyle=' ')
                axes[i].set_ylabel(df_values.columns[i])

                #Calculate a weighted average terminal value for each item
                wavg, wavg_err = weightedAverage(df_values[df_values.columns[i]][-7:-2], df_sigmas[df_sigmas.columns[i]][-7:-2])

                print('Average terminal value for {} = {} +- {}'.format(df_values.columns[i], wavg, wavg_err))
                
            axes[0].set_title(plotTitle)
            axes[-1].set_xlabel('Nucleocentric Distance (arcsec)')


        
        else:
            axes.errorbar(df['Distance (arcsec)'][:-1],df_values[df_values.columns[0]][:-1],yerr=df_sigmas[df_sigmas.columns[0]][:-1],color='C0',marker='x',capsize=2,linestyle=' ')
            axes.set_ylabel(df_values.columns[0])
            axes.set_title(plotTitle)
            axes.set_xlabel('Nucleocentric Distance (arcsec)')

                        #Calculate a weighted average terminal value for each item
            wavg, wavg_err = weightedAverage(df_values[df_values.columns[0]][-7:-2], df_sigmas[df_sigmas.columns[0]][-7:-2])

            print('Average terminal value for {} = {} +- {}'.format(df_values.columns[0], wavg, wavg_err))

        plt.savefig(plotFile)