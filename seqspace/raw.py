# Mapping object for holding raw data for transformed maps
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

from seqspace.base import BaseMap
from seqspace.errors import ErrorMap

class RawMap:
    
    @property
    def phenotypes(self):
        """"""
        return self._phenotypes
        
    @property
    def errors(self):
        """"""
        return self._errors
        
    @phenotypes.setter
    def phenotypes(self, phenotypes):
        """"""
        self._phenotypes = phenotypes
        
    @errors.setter
    def errors(self, errors):
        """"""
        self._errors = ErrorMap()
        self._errors.upper = errors
        self._errors.lower = errors