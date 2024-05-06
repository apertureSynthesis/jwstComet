import numpy as np

def weightedAverage(values,sigmas):
    """
    Calculate a weighted average given a set of values and their uncertainties
    """
    weights = [1/(i**2) for i in sigmas]
    wavg_numerator = np.sum([i*j for i,j in zip(values,weights)])
    wavg_denominator = np.sum(weights)
    wavg_sigma = 1./np.sqrt(np.sum(weights))

    return (wavg_numerator/wavg_denominator), wavg_sigma