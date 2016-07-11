# Mapping object for holding raw data for transformed maps
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

from seqspace.base import BaseMap

class RawMap:
    """ Object for containing untransformed data.
    """
    def __init__(self, GPM):
        self._GPM = GPM

    @property
    def genotypes(self):
        """"""
        return self._genotypes

    @property
    def phenotypes(self):
        """"""
        return self._phenotypes

    @property
    def stdeviations(self):
        """"""
        return self._stdeviations

    @property
    def errors(self):
        """"""
        return self._errors

    @genotypes.setter
    def genotypes(self, genotypes):
        """"""
        self._genotypes = genotypes

    @phenotypes.setter
    def phenotypes(self, phenotypes):
        """"""
        self._phenotypes = phenotypes

    @stdeviations.setter
    def stdeviations(self, stdeviations):
        """"""
        self._stdeviations = stdeviations
