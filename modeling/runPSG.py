import os,sys
import numpy as np
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
from jwstComet.utils import cfgHelper
from jwstComet.modeling import plotPSG

class runPSG(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def getModels(self, specFile, resFile, name, objectType, composition=None, retrieval=None, mode=None, withPlots=False, key=None):

        wave  = []
        spec  = []
        noise = []

        #Read in header parameters, then data
        with open(specFile, 'r') as fn:
            for _, line in enumerate(fn):
                if '#Obs. Start' in line:
                    dateBeg = line.split()[-2] + ' ' + line.split()[-1]
                if '#Obs. End' in line:
                    dateEnd = line.split()[-2] + ' ' + line.split()[-1]

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

        cfgHelper.ephCFG(specFile,name,objectType,midtime,delta,key)

        if retrieval != None:
            cfgHelper.atmCFG(specFile,resFile,composition,retrieval,mode,key)

        if withPlots:
            plotPSG.makePlots(resFile,withPlots)
     



