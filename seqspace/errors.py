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

    def __init__(self, phenotypes, stdeviations, log_transform=False, logbase=np.log10):
        """ If a lower bound is given, use it instead of -variances. """

        self.phenotypes = phenotypes
        self.stdeviations = stdeviations
        self.log_transform = log_transform
        self.logbase = logbase

    @staticmethod
    def transform_upper(bounds, phenotypes, logbase):
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
        return abs(logbase((phenotypes + bounds) /phenotypes))

    @staticmethod
    def transform_lower(bounds, phenotypes, logbase):
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
        return abs(logbase( phenotypes / (phenotypes - bounds) ))


    def wrapper(self, bound, **kwargs):
        """ Wrapper function that changes variances to whatever bound desired. """
        raise Exception(""" Must be implemented in a subclass """)

    @property
    def upper(self):
        """"""
        if self.log_transform:
            return self.transform_upper(self.wrapper(self.stdeviations), self.phenotypes, logbase=self.logbase)
        else:
            return self.wrapper(self.stdeviations)

    @property
    def lower(self):
        """"""
        if self.log_transform:
            return self.transform_lower(self.wrapper(self.stdeviations), self.phenotypes, logbase=self.logbase)
        else:
            return self.wrapper(self.stdeviations)


class StandardDeviationMap(BaseErrorMap):

    def __init__(self, phenotypes, stdeviations, log_transform=False, logbase=np.log10):
        """ Initialize a standard deviations map. """
        super(StandardDeviationMap, self).__init__(phenotypes, stdeviations, log_transform=log_transform, logbase=logbase)

    def wrapper(self, bounds, **kwargs):
        """ Wrapper function to convert Variances if necessary"""
        return bounds

class StandardErrorMap(BaseErrorMap):

    def __init__(self, phenotypes, stdeviations, log_transform=False, n_replicates=2, logbase=np.log10):
        """ Initialize a standard error map object """
        super(StandardErrorMap, self).__init__(phenotypes, stdeviations, log_transform=log_transform, logbase=logbase)
        self.n_replicates = n_replicates

    def wrapper(self, bounds):
        """ Wrapper function to convert Variances if necessary"""
        return bounds/np.sqrt(self.n_replicates)
