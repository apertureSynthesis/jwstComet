import os,sys
from astroquery.jplhorizons import Horizons
from datetime import datetime, timedelta

class cfgBase(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def getCFG(specFile, name, type, mode):

        wave  = []
        spec  = []
        noise = []

        #Read in header parameters, then data
        with open(specFile) as fn:
            for _, line in enumerate(fn):
                if '#Obs. Start' in line:
                    dateBeg = line
                if '#Obs. End' in line:
                    dateEnd = line
                if '#Lower wavelength (um)' in line:
                    waveLo = line.split(':')[1]
                if '#Upper wavelength (um)' in line:
                    waveUp = line.split(':')[1]
                if '#Inner annulus radius (arcsec)' in line:
                    innerRadius = line.split(':')[1]
                if '#Outer annulus radius (arcsec)' in line:
                    outerRadius = line.split(':')[1]
                if '#Aperture radius (arcsec)' in line:
                    radAp = line.split(':')[1]
                if '#X offset (arcsec)' in line:
                    xOffset = line.split(':')[1]
                if '#Y offset (arcsec)' in line:
                    yOffset = line.split(':')[1]
                if '#Wave (micron) Flux (Jy) Noise (Jy)' in line:
                    break
            
            for _, line in enumerate(fn):
                w, s, n = line.split()
                wave.append(w)
                spec.append(s)
                noise.append(n)

        #Run an ephemeris and get pertinent information
        #Find the midpoint of the observations
        start = datetime.strptime(dateBeg,format='%Y-%m-%d %H%M%S')


