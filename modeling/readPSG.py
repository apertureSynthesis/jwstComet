class readPSG(object):
    """
    Read a PSG returned retrieval file

    Inputs
        resFile - ASCII file containing the PSG retrieval results

    Outputs
        returns list of retrieved variables, values, and uncertainties contained in the file
    """


    def __init__(self, resFile):
        super().__init__()
        self.name = self.__class__.__name__
        self.resFile = resFile
        self.retrieval_variables, self.retrieval_values, self.retrieval_sigmas = self.readResults()

    def readResults(self):
        with open(self.resFile, 'r') as fn:
            for index, line in enumerate(fn):
                #Read through and search for the retrieved variables, values, and sigmas
                if '<RETRIEVAL-VARIABLES>' in line:
                    ret_vars = line.split('>')
                    retrieval_variables = ret_vars[-1].split(',')
                if '<RETRIEVAL-VALUES>' in line:
                    ret_vals = line.split('>')
                    retrieval_values = ret_vals[-1].split(',')
                if '<RETRIEVAL-SIGMAS>' in line:
                    ret_sigs = line.split('>')
                    retrieval_sigmas = ret_sigs[-1].split(',')

        for var, val, sig in zip(retrieval_variables, retrieval_values, retrieval_sigmas):
            print('{} = {} +- {}'.format(var,val,sig))

        return [s.strip() for s in retrieval_variables], [float(s.strip()) for s in retrieval_values], [float(s.strip()) for s in retrieval_sigmas]       