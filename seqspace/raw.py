# Mapping object for holding raw data for transformed maps
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

from seqspace.base import BaseMap
from seqspace.errors import StandardErrorMap, StandardDeviationMap

class RawMap(BaseMap):
    """ Object for containing untransformed data.

    Parameters
    ----------
    GPM : GenotypePhenotypeMap object
        The genotype phenotype map object.

    Attributes
    ----------
    genotypes : np.array
        array of binary genotype strings, ordered the same as input from GPM.
    phenotypes : np.array
        phenotypes given by GPM, in the same order as GPM.
    n_replicates : int
        number of replicates.
    logbase : callable
        function for log transforming an array or value.
    stdeviations : array
        standard deviations of genotype phenotype map.
    """
    def __init__(self, GPM):
        self._GPM = GPM
        self.log_transform = False
        self.err = StandardErrorMap(self)
        self.std = StandardDeviationMap(self)

    # -------------------------------------------
    # Getter methods
    # -------------------------------------------

    @property
    def n_replicates(self):
        """Get number of replicates"""
        return self._GPM.n_replicates

    @property
    def logbase(self):
        """Get logbase."""
        return self._GPM.logbase

    @property
    def stdeviations(self):
        """"""
        return self._GPM.stdeviations

    @property
    def genotypes(self):
        """Get genotypes."""
        return self._GPM.genotypes

    @property
    def phenotypes(self):
        """Get non-log-transformed phenotypes"""
        return self._phenotypes

    # -------------------------------------------
    # Setter methods
    # -------------------------------------------

    @phenotypes.setter
    def phenotypes(self, phenotypes):
        """Set the raw phenotypes."""
        self._phenotypes = phenotypes
