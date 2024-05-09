from astropy.io import fits

class readHeader(object):
    """
    Read a JWST header to parse observing information for later storage/use.
    """

    def __init__(self, cubeFile):
        self.cubeFile = cubeFile
        self.target, self.dateBeg, self.dateEnd, self.timeBeg, self.timeEnd, self.setting = self.returnInfo()

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

            setting = detector+'/'+channel+'/'+band

        return target, dateBeg, dateEnd, timeBeg, timeEnd, setting





        