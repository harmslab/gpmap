import numpy as np
import itertools as it

from gpmap.gpm import GenotypePhenotypeMap
from gpmap import utils
from .base import random_mutation_set, BaseSimulation


class NKSimulation(BaseSimulation):
    """Generate genotype-phenotype map from NK fitness model. Creates a table
    with binary sub-sequences that determine the order of epistasis in the
    model.

    The NK fitness landscape is created using a table with binary, length-K,
    sub-sequences mapped to random values. All genotypes are binary with
    length N. The fitness of a genotype is constructed by summing the values
    of all sub-sequences that make up the genotype using a sliding window
    across the full genotype.

    For example, imagine an NK simulation with N=5 and K=2. To construct
    the fitness for the 01011 genotype, select the following sub-sequences
    from an NK table "01", "10", "01", "11", "10". Sum their values together.

    Parameters
    ----------

    Attributes
    ----------
    nk_table : dict
        table with binary sub-sequences as keys which are used to construct
        phenotypes following an NK routine
    keys : array
        array of keys in NK table.
    values : array
        array of values in the NK table.
    """

    def __init__(self, wildtype, mutations, K, k_range=(0, 1),
                 *args, **kwargs):
        super(NKSimulation, self).__init__(wildtype, mutations,
                                           *args, **kwargs)
        # Set parameters
        self.set_order(K)
        self.set_random_values(k_range=k_range)
        self.build()

    @property
    def nk_table(self):
        """NK table mapping binary sequence to value."""
        return self.map("keys", "values")

    @property
    def keys(self):
        """NK table keys.
        """
        return self._keys

    @property
    def values(self):
        """NK table values
        """
        return self._values

    def set_order(self, K):
        """Set the order (K) of the NK model.
        """
        self.K = K
        # point to order
        self.order = self.K
        self._keys = np.array(["".join(r) for r in
                               it.product('01', repeat=self.K)])
        # Reset phenotypes
        self.data['phenotypes'] = np.empty(self.n, dtype=float)

    def set_random_values(self, k_range=(0, 1)):
        """Set the values of the NK table by drawing from a uniform
        distribution between the given k_range.
        """
        if hasattr(self, "keys") is False:
            raise Exception("Need to set K first. Try `set_order` method.")
        self._values = np.random.uniform(k_range[0], k_range[1],
                                         size=len(self.keys))
        self.build()

    def set_table_values(self, values):
        """Set the values of the NK table from a list/array of values.
        """
        if len(values) != len(self.keys):
            raise Exception("Length of the values do not equal the length of "
                            "NK keys. "
                            "Length of keys is : %d" % (len(self.keys),))
        self._values = values
        self.build()

    def build(self):
        """Build phenotypes from NK table.
        """
        nk_table = self.nk_table
        # Check for even interaction
        neighbor = int(self.order / 2)
        if self.order % 2 == 0:
            pre_neighbor = neighbor - 1
        else:
            pre_neighbor = neighbor
        # Use NK table to build phenotypes
        phenotypes = np.zeros(self.n, dtype=float)
        for i in range(len(self.genotypes)):
            f_total = 0
            for j in range(self.length):
                if j - pre_neighbor < 0:
                    pre = self.binary[i][-pre_neighbor:]
                    post = self.binary[i][j:neighbor + j + 1]
                    f = "".join(pre) + "".join(post)
                elif j + neighbor > self.length - 1:
                    pre = self.binary[i][j - pre_neighbor:j + 1]
                    post = self.binary[i][0:neighbor]
                    f = "".join(pre) + "".join(post)
                else:
                    f = "".join(
                        self.binary[i][j - pre_neighbor:j + neighbor + 1])
                f_total += nk_table[f]
            phenotypes[i] = f_total
        self.data.phenotypes = phenotypes
