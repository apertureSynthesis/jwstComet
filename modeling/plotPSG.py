import matplotlib.pyplot as plt
import numpy as np

def makePlots(retFile,withPlots):
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

        #Set up the plots
        fig, axes = plt.subplots(3,1,figsize=(10,10))
        fig.subplots_adjust(hspace=0.25,wspace=0.01)
        axes[0].plot(wave,spec,label='Data')
        axes[0].plot(wave,model,label='Model')
        axes[0].plot(wave,base,label='Base')
        axes[0].legend()

        axes[1].plot(wave,spec-model,label='Residual')
        axes[1].plot(wave,base,label='Base')
        axes[1].legend()

        axes[2].plot(wave,model-base,label='Model Gas Spectra')
        axes[2].legend()
        plt.show()