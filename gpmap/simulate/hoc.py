from .nk import NKSimulation
from .base import random_mutation_set


class HouseOfCardsSimulation(NKSimulation):
    """Construct a 'House of Cards' fitness landscape.
    """

    def __init__(self, wildtype, mutations, k_range=(0, 1), *args, **kwargs):
        super(NKSimulation, self).__init__(wildtype, mutations, *args,
                                           **kwargs)
        # Set parameters
        self.set_order(len(self.binary[0]))
        self.set_random_values(k_range=k_range)
        self.build()
