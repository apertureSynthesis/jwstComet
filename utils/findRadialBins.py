import numpy as np
from astropy.io import fits
import astropy.units as u
from jwstComet.utils import readCube

class findRadialBins(object):

    def __init__(self, cubefile, annulusWidth=0.1*u.arcsec):
        self.cubefile = cubefile
        self.annulusWidth = annulusWidth
        self.r_ins, self.r_outs, self.r_mids = self.generateAnnuli()

    def generateAnnuli(self):
        
        #Get the geometry of the cube
        sciCube = readCube(self.cubefile)

        distances = []

        #Find the most distant pixel on the chip from the photocenter
        for x in range(sciCube.xs):
            for y in range(sciCube.ys):
                #Check that we're on the chip
                if sciCube.wmap[50,y,x] != 0:
                    #Calculate distance
                    dx = x - sciCube.xcenter
                    dy = y - sciCube.ycenter
                    dr = np.sqrt(dx**2 + dy**2)
                    distances.append(dr)

        #Find the largest distance, only go about 80% of the way
        #to make sure we mostly remain on the chip
        rad_max = np.round(0.8*np.nanmax(distances))
        

        #Work out how big our annulus width is in pixels
        pixar_a2 = sciCube.hdr['PIXAR_A2']*u.arcsec**2
        pixscale = u.pixel_scale(np.sqrt(pixar_a2) / u.pixel)
        annulusWidthPix = (self.annulusWidth).to(u.pixel,pixscale)
        
        #Generate the inner and outer radii
        r_ins = np.arange(0, rad_max, annulusWidthPix.value)
        r_outs = np.arange(annulusWidthPix.value, rad_max + annulusWidthPix.value, annulusWidthPix.value)
        r_mids = [(i+j)/2. for i,j in zip(r_outs,r_ins)]

        #Convert radii to arcseconds
        r_ins_as  = [(r_in*u.pixel).to(u.arcsec,pixscale) for r_in in r_ins]
        r_outs_as = [(r_out*u.pixel).to(u.arcsec,pixscale) for r_out in r_outs]
        r_mids_as = [(r_mid*u.pixel).to(u.arcsec,pixscale) for r_mid in r_mids]

        return r_ins_as, r_outs_as, r_mids_as


