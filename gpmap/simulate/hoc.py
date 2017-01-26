from .nk import NKSimulation

class HouseOfCardsSimulation(NKSimulation):
    """Construct a 'House of Cards' fitness landscape.
    """
    @classmethod
    def from_length(cls, length, k_range=(0,1), *args, **kwargs):
        """Build an House of cards model form

        Parameters
        ----------
        length : int
            length of genotypes
        K : int
            Order of epistasis in NK model.

        Returns
        -------
        self : NKSimulation
        """
        mutations = dict([(i,["0","1"]) for i in range(length)])
        wildtype = "".join([m[0] for m in mutations.values()])
        self = cls(wildtype, mutations, length, k_range=k_range, *args, **kwargs)
        return self
