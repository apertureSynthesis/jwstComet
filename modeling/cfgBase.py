import os,sys
from astroquery.jplhorizons import Horizons
from jwstComet.utils import cfgHelper

class cfgBase(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def getCFG(self,specFile, name, objectType, composition=None, retrieval=None, mode=None):

        wave  = []
        spec  = []
        noise = []

        #Read in header parameters, then data
        with open(specFile) as fn:
            for _, line in enumerate(fn):
                if '#Obs. Start' in line:
                    dateBeg = line.split()[-2] + ' ' + line.split()[-1]
                if '#Obs. End' in line:
                    dateEnd = line.split()[-2] + ' ' + line.split()[-1]
                if '#Lower wavelength (um)' in line:
                    waveLo = line.split()[-1]
                if '#Upper wavelength (um)' in line:
                    waveUp = line.split()[-1]
                if '#Inner annulus radius (arcsec)' in line:
                    innerRadius = line.split()[-1]
                if '#Outer annulus radius (arcsec)' in line:
                    outerRadius = line.split()[-1]
                if '#Aperture radius (arcsec)' in line:
                    radAp = line.split()[-1]
                if '#X offset (arcsec)' in line:
                    xOffset = line.split()[-1]
                if '#Y offset (arcsec)' in line:
                    yOffset = line.split()[-1]
                if '#Wave (micron) Flux (Jy) Noise (Jy)' in line:
                    break
            
            for _, line in enumerate(fn):
                w, s, n = line.split()
                wave.append(w)
                spec.append(s)
                noise.append(n)

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



