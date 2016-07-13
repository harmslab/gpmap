# Mapping Object for mutations-to-genotypes for genotype-phenotype maps
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

from seqspace.base import BaseMap

class MutationMap(BaseMap):
    """This object tracks the index and order of mutations in an epistatic map.

    Parameters
    ----------
    GPM : GenotypePhenotypeMap object
        The genotype phenotype map object to translate as Binary.

    Attributes
    ----------
    wildtype : str
        wildtype genotype
    mutations : dict
        mapping of site number to alphabet
    n : int
        number of mutations
    """
    def __init__(self, GPM):
        self._GPM = GPM

    # ------------------------------------------------------------------
    # Getter methods for attributes that are not set explicitly by user.
    # ------------------------------------------------------------------

    @property
    def wildtype(self):
        """ Get possible that occur from reference system. """
        return self._GPM.wildtype

    @property
    def mutant(self):
        """Get the farthest mutant in genotype-phenotype map."""
        _mutant = []
        _wt = self.wildtype
        for i in range(0,self.n):
            site = _wt[i]
            options = self.mutations[i]
            for o in options:
                if o != site:
                    _mutant.append(o)
        return "".join(_mutant)

    @property
    def mutations(self):
        """ Get possible that occur from reference system. """
        return self._mutations

    @property
    def n(self):
        """ Get the number of mutations in the space. """
        return self._n

    # ------------------------------------------------------------------
    # Setter methods for attributes that are not set explicitly by user.
    # ------------------------------------------------------------------

    @mutations.setter
    def mutations(self, mutations):
        """ Set the mutation alphabet for all sites in wildtype genotype.

        Examples
        --------
        mutations = { site_number : alphabet }. If the site
        alphabet is note included, the model will assume binary
        between wildtype and derived::

            mutations = {
                0: [alphabet],
                1: [alphabet],

            }
        """
        if type(mutations) != dict:
            raise TypeError("mutations must be a dict")
        self._mutations = mutations
        self._n = len(mutations)
