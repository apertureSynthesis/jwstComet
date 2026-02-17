import matplotlib.pyplot as plt
import numpy as np

def makePlots(retFile,withPlots=False):
    """
    Reads a PSG retrieval result file. Provides a model-data comparison plot if desired

    Inputs
        retFile - ASCII file containing the retrieval results
        withPlots - whether to plot the model-data comparison

    Outputs
        If desired, shows (but not saves) plots
    """
    #Plot the model-data fit if desired
    if withPlots:
        wave  = []
        spec  = []
        dspec = []
        model = []
        base  = []
        #Read the file line by line
        with open(retFile, 'r') as fn:
            for line in fn:
                if 'results_dat.txt' in line:
                    break

            for line in fn:
                if 'results_log.txt' in line:
                    break
                hasBase = True
                try:
                    lwave,lspec,ldspec,lmodel,lbase = line.split()
                    wave = np.concatenate([wave,[float(lwave)]])
                    spec = np.concatenate([spec,[float(lspec)]])
                    dspec = np.concatenate([dspec,[float(ldspec)]])
                    model = np.concatenate([model,[float(lmodel)]])
                    base = np.concatenate([base,[float(lbase)]])
                except:
                    hasBase = False
                    lwave,lspec,ldspec,lmodel = line.split()
                    wave = np.concatenate([wave,[float(lwave)]])
                    spec = np.concatenate([spec,[float(lspec)]])
                    dspec = np.concatenate([dspec,[float(ldspec)]])
                    model = np.concatenate([model,[float(lmodel)]])                   

        #Set up the plots
        #try:
        fig, axes = plt.subplots(4,1,figsize=(10,10))
        fig.subplots_adjust(hspace=0.25,wspace=0.01)
        if hasBase:
            print(f'Band sum = {np.nansum(model-base):.3e}')
        axes[0].plot(wave,spec,label='Data')
        axes[0].plot(wave,model,label='Model')
        #axes[0].plot(wave,-1*base,label='Base')
        axes[0].legend()

        axes[1].plot(wave,spec-model,label='Residual')
        #axes[1].plot(wave,-1*base,label='Base')
        axes[1].legend()

        if hasBase:
            axes[2].plot(wave,model-base,label='Model Gas Spectra')
        else:
            axes[2].plot(wave,model,label='Model Gas Spectra')
        axes[2].legend()

        if hasBase:
            axes[3].plot(wave,-1*base,label='Base')
        else:
            axes[3].axhline(0,label='Base')
        axes[3].legend()
        plt.show()
       # except:
        #    print('Missing spectral fit')
       #     pass