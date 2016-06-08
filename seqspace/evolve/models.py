import numpy as np

def fixation(fitness1, fitness2, N=10e8):
    """ Simple fixation probability between two organism with fitnesses 1 and 2.

    Note that N is the effective population size.
    """
    sij = (fitness2 - fitness1)/abs(fitness1)

    # Check the value of denominator
    denominator = 1 - np.exp(-N * sij)
    numerator = 1 - np.exp(- sij)

    # Calculate the fixation probability
    fixation = numerator / denominator
    return fixation
