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
        AA = utils.AMINO_ACIDS
        random.shuffle(AA)
        alphabet = AA[:size[i]]
        mutations[i] = alphabet
    return mutations


class BaseSimulation(GenotypePhenotypeMap):
    """ Build a simulated GenotypePhenotypeMap. Generates random phenotypes.
    """

    def __init__(self, wildtype, mutations, *args, **kwargs):
        # build genotypes
        genotypes = utils.mutations_to_genotypes(wildtype, mutations)
        phenotypes = np.empty(len(genotypes), dtype=float)
        super(BaseSimulation, self).__init__(wildtype,
                                             genotypes,
                                             phenotypes,
                                             mutations=mutations,
                                             *args,
                                             **kwargs)

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

    def set_stdeviations(self, sigma):
        """Add standard deviations to the simulated phenotypes, which can then be
        used for sampling error in the genotype-phenotype map.

        Parameters
        ----------
        sigma : float or array-like
            Adds standard deviations to the phenotypes. If float, all
            phenotypes are given the same stdeviations. Else, array must be
            same length as phenotypes and will be assigned to each phenotype.
        """
        stdeviations = np.ones(len(self.phenotypes)) * sigma
        self.data.stdeviations = stdeviations
        return self

    def build(self):
        raise Exception("must be implemented in subclass.")
