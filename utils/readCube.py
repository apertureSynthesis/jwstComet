import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip

class readCube(object):
    """
    Reads a JWST IFU cube. Returns key quantities from the cube, including data and header information.

    Inputs
        cubeFile - FITS file containing the datacube

    Outputs
        data - 3D spectral cube data
        derr - 3D uncertainty cube
        qty - 3D quality cube
        wmap - 3D cube for whether pixel is on the detector
        cdata - 2D spectrally integrated cube data
        xcenter - horizontal pixel coordinate of photocenter
        xs - horizontal size of the array
        ycenter - vertical pixel coordinate of the photocenter
        ys - vertical size of the array
        hdr - header information from the cube

    """

    def __init__(self, cubeFile):
        super().__init__()
        self.name = self.__class__.__name__
        self.cubeFile = cubeFile
        self.data, self.derr, self.qty, self.wmap, self.cdata, self.xcenter, self.xs, self.ycenter, self.ys, self.hdr = self.returnCube()

    def returnCube(self):
        #Open the file and read the data
        with fits.open(self.cubeFile) as fr:
            hdr = fr["sci"].header
            data = fr["sci"].data
            derr = fr["err"].data
            qty  = fr["dq"].data
            wmap = fr["wmap"].data

        xs = data.shape[2]
        ys = data.shape[1]

        #Mask low-quality points
        mm = np.nanmedian(derr)
        upper_threshold = 1e4*mm
        lower_threshold = -1e3*mm
        ind = (qty!=0) + (data>upper_threshold) + (data<lower_threshold)
        data[ind] = np.nan

        #Collapse the cube and find the photocenter
        cdata = np.nansum(data, axis=0)
        #Set pixels off chip to nan
        k = np.where(wmap[50,:,:] == 0)
        cdata[k] = np.nan
        ind_max = np.where(cdata == np.nanmax(cdata))
        xcenter = ind_max[1][0]
        ycenter = ind_max[0][0]  

        return data, derr, qty, wmap, cdata, xcenter, xs, ycenter, ys, hdr     