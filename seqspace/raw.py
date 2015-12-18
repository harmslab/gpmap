# Mapping object for holding raw data for transformed maps
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

from seqspace.base import BaseMap

class RawMap:
    
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
        
    @phenotypes.setter
    def phenotypes(self, phenotypes):
        """"""
        self._phenotypes = phenotypes
        
    @stdeviations.setter
    def stdeviations(self, stdeviations):
        """"""
        self._stdeviations = stdeviations
