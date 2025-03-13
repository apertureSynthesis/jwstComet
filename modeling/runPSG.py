import os,sys
import numpy as np
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
from jwstComet.utils import cfgHelper
from jwstComet.modeling import plotPSG

class runPSG(object):

    """
    Lay out the process for running the requested PSG models
    """

    def __init__(self):
        super().__init__()
        self.name = self.__class__.__name__

    def getModels(self, specFile, resFile, name, objectType, composition=None, retrieval=None, mode=None, withCont=False, withPlots=False, withEph=True, local=True):
        
        """
        Run a requested PSG model and retrieval.

        Inputs
            specFile - ASCII file containing the spectrum to be analyzed
            resFile - ASCII file for saving the returned results from the PSG
            name - name of the comet or asteroid
            objectType - type of small body: asteroid or comet
            composition - dictionary containing compositional information for building the PSG model atmosphere. optional
            retrieval - dictionary containing quantities to be retrieved for each PSG model run. optional.
            mode - extraction mode (circle, rectangle, mapping, azimuthal)
            withCont - whether we are asking the PSG to simulate the continuum or instead simply subtract a baseline
            withPlots - whether we are plotting the results
            withEph - whether we are asking the PSG to retrieve ephemeris parameters or are instead using a local copy
            local - are we interrogating a local copy of the PSG or instead sending requests to the online server
        """


        #Read in header parameters for the observation times
        with open(specFile, 'r') as fn:
            for _, line in enumerate(fn):
                if '#Obs. Start' in line:
                    dateBeg = line.split()[-2] + ' ' + line.split()[-1]
                if '#Obs. End' in line:
                    dateEnd = line.split()[-2] + ' ' + line.split()[-1]

        if withEph:
            #Run an ephemeris and get pertinent information
            #Find the midpoint of the observations
            obj = Horizons(id = name, id_type = 'designation', location='@JWST',
                        epochs = {'start':dateBeg, 'stop':dateEnd, 'step':'1m'})
            eph = obj.ephemerides(quantities = '19,20,23,24', no_fragments=True, closest_apparition=True)

            df_eph = eph.to_pandas()

            obs_midpoint = int(len(df_eph)/2.)
            delta = df_eph['delta'][obs_midpoint]
            midtime = df_eph['datetime_str'][obs_midpoint]
            print(midtime)

            cfgHelper.ephCFG(specFile,name,objectType,midtime,local)

        #Run a retrieval if requested
        if retrieval != None:
            cfgHelper.atmCFG(specFile,resFile,composition,retrieval,mode,withCont,local)

        if withPlots:
            plotPSG.makePlots(resFile,withPlots)
     



