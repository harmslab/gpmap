import numpy as np
from .base import BaseSimulation


class RandomPhenotypesSimulation(BaseSimulation):
    """ Build a simulated GenotypePhenotypeMap. Generates random phenotypes.
    """

    def __init__(self, wildtype, mutations, phenotype_range=(0, 1),
                 *args, **kwargs):
        # build genotypes
        super(RandomPhenotypesSimulation, self).__init__(wildtype,
                                                         mutations,
                                                         *args, **kwargs)
        self.phenotype_range = phenotype_range
        self.build()

    def build(self):
        """Build phenotypes"""
        low, high = self.phenotype_range[0], self.phenotype_range[1]
        self.data['phenotypes'] = np.random.uniform(low, high,
                                                    size=len(self.genotypes))
