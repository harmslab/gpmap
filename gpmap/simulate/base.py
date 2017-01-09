import random
from gpmap import utils

def random_mutation_set(length, alphabet_size=2):
    """Generate a random mutations dictionary for simulations.

    Parameters
    ----------
    length : length of genotypes

    alphabet_size : int or list
        alphabet size at each site. if list is given, will make site i have
        size alphab_size[i].
    """
    if type(alphabet_size) == int:
        size = [alphabet_size for i in range(length)]
    else:
        size = alphabet_size

    alphabet = utils.AMINO_ACIDS[:size]
    mutations = dict([(i, alphabet[i])for i in range(length)])
    return mutations
