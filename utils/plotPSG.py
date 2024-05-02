import matplotlib.pyplot as plt
import numpy as np

def plotPSG(retFile,withPlots):
    """
    Reads a PSG retrieval result file. Provides a model-data comparison plot if desired
    """
    #Plot the model-data fit if desired
    if withPlots:
        wave  = []
        spec  = []
        dspec = []
        model = []
        base  = []
        with open(retFile, 'r') as fn:
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