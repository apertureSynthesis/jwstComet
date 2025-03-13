from astropy.io import fits

class readHeader(object):
    """
    Read a JWST header to parse observing information for later storage/use.

    Inputs
        cubeFile - cubeFile - FITS file containing the datacube

    Outputs - information to be stored in the header of extracted spectra files
        target - target name 
        dateBeg - UT beginning date of observations 
        dateEnd - UT end date of observations 
        timeBeg - UT beginning time of observations 
        timeEnd - UT end time of observations 
        instrument - name of JWST instrument used (only MIRI, NIRSpec currently supported) 
        setting - name of instrumental setting used

    """

    def __init__(self, cubeFile):
        super().__init__()
        self.name = self.__class__.__name__
        self.cubeFile = cubeFile
        self.target, self.dateBeg, self.dateEnd, self.timeBeg, self.timeEnd, self.instrument, self.setting = self.returnInfo()

    def returnInfo(self):
        #Open the file and read the Primary header
        with fits.open(self.cubeFile) as fr:
            info = fr[0].header

        #target name
        target = info['TARGNAME']

        #beginning/end datetime
        timeStampBeg = info['DATE-BEG']
        timeStampEnd = info['DATE-END']
        
        dateBeg = timeStampBeg.split('T')[0]
        timeBeg = timeStampBeg.split('T')[1]

        dateEnd = timeStampEnd.split('T')[0]
        timeEnd = timeStampEnd.split('T')[1]

        instrument = info['INSTRUME']
        if instrument == 'NIRSPEC':

            grating = info['GRATING']
            filter = info['FILTER']

            setting = grating+'/'+filter

        if instrument == 'MIRI':

            detector = info['DETECTOR']
            channel = info['CHANNEL']
            band = info['BAND']

            setting = channel+'/'+band.strip()

        return target, dateBeg, dateEnd, timeBeg, timeEnd, instrument, setting





        