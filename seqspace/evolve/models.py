import numpy as np

def fixation(fitness1, fitness2):
    """ Simple fixation probability between two organism with fitnesses 1 and 2. """
    sij = (fitness2 - fitness1)/abs(fitness1)
    
    # Set negative probabilities to 0
    if sij < 0:
        sij = 0
        
    fixation = 1 - np.exp(-2*sij)
    return fixation
