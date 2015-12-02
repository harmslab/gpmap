# Mapping object for holding upper and lower error bars
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

from seqspace.base import BaseMap

class ErrorMap:
    
    @property
    def upper(self):
        """"""
        return self._upper
        
    @property
    def lower(self):
        return self._lower
        
    @upper.setter
    def upper(self, upper):
        """"""
        self._upper = upper
        
    @lower.setter
    def lower(self, lower):
        """"""
        self._lower = lower