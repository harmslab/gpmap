import numpy as np

def fixation(fitness1, fitness2, N=10e8, *args, **kwargs):
    """ Simple fixation probability between two organism with fitnesses 1 and 2.

    Note that N is the effective population size.

    .. math::
        p_{\\text{fixation}} = \\frac{1 - e^{-N \\frac{f_2-f_1}{f1}}}{1 - e^{-\\frac{f_2-f_1}{f1}}}
    """
    sij = (fitness2 - fitness1)/abs(fitness1)
    # Check the value of denominator
    denominator = 1 - np.exp(-N * sij)
    numerator = 1 - np.exp(- sij)
    # Calculate the fixation probability
    if denominator == 0:
        fixation = 0
    else:
        fixation = numerator / denominator
    return fixation

def wright_fisher(fitness1, fitness2, N=10e8, *args, **kwargs):
    """ Calculates a fixation probability using the Wright-Fisher process
    (a birth-death process).

    .. math::
        p_{\\text{fixation}} = \\frac{1 - (\\frac{f_1}{f2})^{2}}{1 - (\\frac{f_1}{f2})^{2N}}
    """
    if fitness2 == 0:
        return 0
    else:
        frac = fitness1 / fitness2
        if frac == 1:
            return 0
        else:
            return (1 - frac ** 2) / (1 - frac**(2*N))
