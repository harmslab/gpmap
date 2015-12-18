# Mapping object for holding upper and lower error bars
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

import numpy as np

from seqspace.base import BaseMap


class BaseErrorMap(BaseMap):
    
    def __init__(self, phenotypes, variances, log_transform=False, lower=None):
        """ If a lower bound is given, use it instead of -variances. """
        
        self.phenotypes = phenotypes
        self.variances = variances
        self.log_transform = log_transform
    
        # Set a lower bound -- if None, just use upper bound.
        self.variances_ = variances
        if lower is None:
            self.variances_ = variances
        
        
    
    @staticmethod
    def transform(bounds, phenotypes):
        """ Log transformation scaling. 
        
            Untransformed data looks as so:
                
                Yupper = Ymean + bound
                Ylower = Ymean - bound
                
            We want log(bounds)
                ie.
                    log(Yupper) - log(Ymean)
                    log(Ylower) + log(Ymean)
                
            so log(bound) = log(1 + bound/Ymean)
               log(bound) = log(1 - bound/Ymean)
        """
        return abs(np.log10(1 + bounds/phenotypes))
        
    @staticmethod
    def wrapper(bound, **kwargs):
        """ Wrapper function that changes variances to whatever bound desired. """
        raise Exception(""" Must be implemented in a subclass """)
    
    @property
    def upper(self):
        """"""
        if self.log_transform:
            return self.transform(self.wrapper(self.variances), self.phenotypes)
        else:
            return self.wrapper(self.variances)
        
    @property
    def lower(self):
        """"""
        if self.log_transform:
            return self.transform(-self.wrapper(self.variances), self.phenotypes)
        else:
            return self.wrapper(self.variances_) 
        
        
class VarianceMap(BaseErrorMap):
    
    def __init__(self, phenotypes, variances, log_transform=False, lower=None):
        """ Initialize a variance map object """
        super(VarianceMap, self).__init__(phenotypes, variances, log_transform=log_transform, lower=None)
    
    @staticmethod
    def wrapper(bounds, **kwargs):
        """ Wrapper function to convert Variances if necessary"""
        return bounds
            

class StandardDeviationMap(BaseErrorMap):
    
    def __init__(self, phenotypes, variances, log_transform=False, lower=None):
        """ Initialize a standard deviations map. """
        super(StandardDeviationMap, self).__init__(phenotypes, variances, log_transform=log_transform, lower=None)

    @staticmethod
    def wrapper(bounds, **kwargs):
        """ Wrapper function to convert Variances if necessary"""
        return np.sqrt(bounds)

class StandardErrorMap(BaseErrorMap):
    
    def __init__(self, phenotypes, variances, log_transform=False, n_replicates=2, lower=None):
        """ Initialize a standard error map object """
        super(StandardErrorMap, self).__init__(phenotypes, variances, log_transform=log_transform, lower=None)
        self.n_replicates = n_replicates
        
    @staticmethod
    def wrapper(bounds, n_replicates=2):
        """ Wrapper function to convert Variances if necessary"""
        return np.sqrt(bounds)/np.sqrt(n_replicates)