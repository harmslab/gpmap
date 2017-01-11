import random
import numpy as np
from gpmap import utils
from gpmap.gpm import GenotypePhenotypeMap

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
    # build mutations dictionary
    mutations = {}
    for i in range(length):
        alphabet = utils.AMINO_ACIDS[:size[i]]
        mutations[i] = alphabet
    return mutations

class GenotypePhenotypeSimulation(GenotypePhenotypeMap):
    """ Build a simulated GenotypePhenotypeMap. Generates random phenotypes.
    """
    def __init__(self, wildtype, mutations, range=(0,1), *args, **kwargs):
        # build genotypes
        genotypes = utils.mutations_to_genotypes(wildtype, mutations)
        phenotypes = np.empty(len(genotypes), dtype=float)
        super(GenotypePhenotypeSimulation, self).__init__(wildtype, genotypes,
            phenotypes,
            *args,
            **kwargs,
        )
        self.set_random(range=range)

    def set_random(self, range=(0,1)):
        """ Get a set of random
        """
        self.phenotypes = np.random.uniform(range[0], range[1], size=self.n)

    @classmethod
    def from_length(cls, length, alphabet_size=2, *args, **kwargs):
        """ Create a simulate genotype-phenotype map from a given genotype length.

        Parameters
        ----------
        length : int
            length of genotypes
        alphabet_size : int (optional)
            alphabet size

        Returns
        -------
        self : GenotypePhenotypeSimulation
        """
        mutations = random_mutation_set(length, alphabet_size=alphabet_size)
        wildtype = "".join([m[0] for m in mutations.values()])
        self = cls(wildtype, mutations, *args, **kwargs)
        return self
