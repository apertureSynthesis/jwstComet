import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
from jwstComet.extraction import beam
from jwstComet.modeling import runPSG
from jwstComet.utils import readCube

class mapping(object):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @u.quantity_input(waveLo=u.um, waveUp=u.um, radAp=u.arcsec)
    def makeMaps(self,cubeFiles,specStem,csvFile,waveLo,waveUp,radAp):
        """
        Read in a JWST IFU cube. Find the photocenter. Extract spectra across the
        entire cube. Send them to the PSG for analysis. Plot the results.
        """

        #Read in the first cube to serve as a coordinate reference
        sciCube = readCube(cubeFiles[0])

        retrieval_x_indexes   = []
        retrieval_x_offsets   = []
        retrieval_y_indexes   = []
        retrieval_y_offsets   = []
        retrieval_variables   = []
        retrieval_values      = []
        retrieval_sigmas      = []
