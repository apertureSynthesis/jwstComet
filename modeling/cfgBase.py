import os,sys
import numpy as np
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
from jwstComet.utils import cfgHelper

class cfgBase(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def getCFG(self, specFile, name, objectType, composition=None, retrieval=None, mode=None, withPlots=False):

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
        obj = Horizons(id = name, id_type = None, location='@JWST',
                       epochs = {'start':dateBeg, 'stop':dateEnd, 'step':'1m'})
        eph = obj.ephemerides(quantities = '19,20,23,24')

        df_eph = eph.to_pandas()

        obs_midpoint = int(len(df_eph)/2.)
        delta = df_eph['delta'][obs_midpoint]
        midtime = df_eph['datetime_str'][obs_midpoint]
        print(midtime)

        cfgHelper.ephCFG(specFile,name,objectType,midtime,delta)
        cfgHelper.atmCFG(specFile,composition,retrieval,mode)

        if withPlots:
            resName = specFile[:-3]+'ret-results.txt'
            wave  = []
            spec  = []
            dspec = []
            model = []
            base  = []
            with open(resName, 'r') as fn:
                for line in fn:
                    if 'results_dat.txt' in line:
                        break

                for line in fn:
                    if 'results_log.txt' in line:
                        break
                    try:
                        lwave,lspec,ldspec,lmodel,lbase = line.split()
                        wave = np.concatenate([wave,[float(lwave)]])
                        spec = np.concatenate([spec,[float(lspec)]])
                        dspec = np.concatenate([dspec,[float(ldspec)]])
                        model = np.concatenate([model,[float(lmodel)]])
                        base = np.concatenate([base,[float(lbase)]])
                    except:
                        pass

            plt.plot(wave,spec,label='Data')
            plt.plot(wave,model,label='Model')
            plt.legend()
            plt.show()            



