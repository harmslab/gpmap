# Mapping object for holding upper and lower error bars
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

import numpy as np

from seqspace.base import BaseMap
#from seqspace.stats import corrected_sterror, corrected_std

class BaseErrorMap(BaseMap):
    
    def __init__(self, upper, lower, func, **kwargs):
        self._upper = upper
        self._lower = lower
        self._kwargs = kwargs
        self._func = func
    
    @property
    def upper(self):
        """"""
        return self._func(self._upper, **self._kwargs)
        
    @property
    def lower(self):
        """"""
        return self._func(self._lower, **self._kwargs)
        
        
class VarianceMap(BaseErrorMap):
    
    def __init__(self, upper, lower):
        """ Initialize a variance map object """
        super(VarianceMap, self).__init__(upper, lower, func=lambda x: x)


class StandardDeviationMap(BaseErrorMap):
    
    def __init__(self, upper, lower, n_replicates=2):
        """ Initialize a standard deviations map. """
        super(StandardDeviationMap, self).__init__(upper, lower, np.sqrt)


class StandardErrorMap(BaseErrorMap):
    
    def __init__(self, upper, lower, n_replicates=2):
        """ Initialize a standard error map object """
        super(StandardErrorMap, self).__init__(upper, lower, lambda x: np.sqrt(x)/np.sqrt(n_replicates))