import numpy as np
import random
from gpmap.gpm import GenotypePhenotypeMap
from gpmap import utils
from .base import random_mutation_set, BaseSimulation



class MultiPeakMountFujiSimulation(BaseSimulation):
    """Constructs a genotype-phenotype map from a Mount Fuji model. [1]_

    A Mount Fuji sets a "global" fitness peak (max) on a single genotype in
    the space. The fitness goes down as a function of hamming distance away
    from this genotype, called a "fitness field". The strength (or scale) of
    this field is linear and depends on the parameters `field_strength`.
    Roughness can be added to the Mount Fuji model using a random
    `roughness` parameter. This assigns a random

    .. math::

        f(g) = \\nu (g) + c \cdot d(g_0, g)

    where $\\nu$ is the roughness parameter, $c$ is the field strength,
    and $d$ is the hamming distance between genotype $g$ and the
    reference genotype.

    Parameters
    ----------
    wildtype : str
        reference genotype to put the

    mutations : dict
        mutations alphabet for each site

    field_strength : float
        field strength

    roughness_width : float
        Width of roughness distribution

    roughness_dist : str, 'normal'
        Distribution used to create noise around phenotypes.


    References
    ----------

    _ [1] Szendro, Ivan G., et al. "Quantitative analyses of empirical fitness
        landscapes." Journal of Statistical Mechanics: Theory and Experiment
        2013.01 (2013): P01005.
    """

    def __init__(
            self,
            wildtype,
            mutations,
            field_strength=1,
            roughness_width=None,
            roughness_dist='normal',
            peak_n=2,
            min_dist=1,
            max_dist=None,
            a_state=None,
            b_state=None,
            *args,
            **kwargs):
        # Call parent class.
        super(MultiPeakMountFujiSimulation, self).__init__(wildtype, mutations,
                                                  *args, **kwargs)
        # Set the field strength and roughness
        self._field_strength = field_strength
        self._roughness_width = roughness_width
        self._roughness_dist = roughness_dist
        self._roughness = None
        self._peak_n = peak_n
        self._min_dist = min_dist
        self._max_dist = max_dist
        self._a_state = a_state
        self._b_state = b_state
        self.build()

    @classmethod
    def from_length(
        cls,
        length,
        field_strength=1,
        roughness_width=None,
        roughness_dist='normal',
        *args,
        **kwargs):
        """Constructs a genotype-phenotype map from a Mount Fuji model. [1]_

        A Mount Fuji sets a "global" fitness peak (max) on a single genotype in
        the space. The fitness goes down as a function of hamming distance away
        from this genotype, called a "fitness field". The strength (or scale) of
        this field is linear and depends on the parameters `field_strength`.
        Roughness can be added to the Mount Fuji model using a random
        `roughness` parameter. This assigns a random

        .. math::

            f(g) = \\nu (g) + c \cdot d(g_0, g)

        where $\\nu$ is the roughness parameter, $c$ is the field strength,
        and $d$ is the hamming distance between genotype $g$ and the
        reference genotype.

        Parameters
        ----------
        length : int
            length of the genotypes.

        field_strength : float
            field strength

        roughness_width : float
            Width of roughness distribution

        roughness_dist : str, 'normal'
            Distribution used to create noise around phenotypes.
        """
        cls = super(MultiPeakMountFujiSimulation, cls).from_length(
            length,
            field_strength=field_strength,
            roughness_width=roughness_width,
            roughness_dist=roughness_dist,
            *args,
            **kwargs
        )
        return cls

    @property
    def a_state(self):
        """Wild type state"""
        if self._a_state is not None:
            return self._a_state
        elif self._a_state is None:
            self._a_state = self.wildtype
            return self._a_state

    @property
    def b_state(self):
        """Derived state"""
        if self._b_state is not None:
            return self._b_state
        elif self._b_state is None:
            self._b_state = utils.farthest_genotype(self.a_state, self.genotypes)
            return self._b_state

    @property
    def peaks(self):
        """Find n peaks that meet the max_dist/min_dist requirement"""
        self._peaks = [self.b_state, self.a_state]
        while len(self._peaks) < self.peak_n:
            proposed = random.choice(self.genotypes)  # Propose a new peak.
            add = False
            for peak in self._peaks:
                if utils.hamming_distance(peak, proposed) >= self.min_dist <= self.max_dist:  # Check dist. requirements
                    add = True
                else:
                    break
            if add:
                self._peaks.append(proposed)
        return self._peaks

    @property
    def hamming(self):
        """Hamming distances from each peak"""
        try:
            return self._hamming
        # calculate the hamming distance if not done already
        except AttributeError:
            hd = np.empty([len(self.peaks), len(self.genotypes)], dtype=int)
            for i, peak in enumerate(self.peaks):
                for j, g in enumerate(self.genotypes):
                    hd[i][j] = utils.hamming_distance(peak, g)
            self._hamming = hd
            return self._hamming

    @property
    def peak_n(self):
        """Number of peaks"""
        return self._peak_n

    @property
    def max_dist(self):
        """Maximum hamming distance between two peaks"""
        if self._max_dist is not None:
            return self._max_dist
        # If maximum distance between peaks is not given, set to the maximum hamming distance.
        elif self._max_dist is None:
            self._max_dist = len(list(self.genotypes[0]))
            return self._max_dist

    @property
    def min_dist(self):
        """Minimum hamming distance between two peaks"""
        return self._min_dist

    @property
    def roughness(self):
        """Array of roughness values for all genotypes"""
        if self._roughness is not None:
            return self._roughness

        elif self.roughness_width is None:
            return np.zeros(self.n)

        elif self.roughness_dist == 'normal':
            # Set roughness.
            self._roughness = np.random.normal(
                scale=self.roughness_width,
                size=self.n)

            return self._roughness

        elif self.roughness_dist == 'uniform':
            # Set roughness.
            self._roughness = np.random.uniform(
                high=self.roughness_width,
                low=-self.roughness_width,
                size=self.n)

            return self._roughness

        else:
            raise Exception("Roughness isn't set.")

    @property
    def roughness_dist(self):
        """Roughness distribution."""
        return self._roughness_dist

    @roughness_dist.setter
    def roughess_dist(self, roughness_dist):
        """Set the roughness distribution. Also sets the roughness array.s
        """
        # Make sure roughness dist is the right type.
        if not isinstance(roughness_dist, str):
            raise TypeError('roughness_dist must be a string.')

        # Get roughness distribution
        if roughness_dist not in ['normal', 'uniform']:
            raise AttributeError('roughness_dist must be '
                                 'either normal or uniform')

        # Set roughness distribution and reset map
        self._roughness_dist = roughness_dist
        self._roughness = None
        self.build()

    @property
    def roughness_width(self):
        return self._roughness_width

    @roughness_width.setter
    def roughness_width(self, roughness_width):
        # Set roughness distribution and reset map
        self._roughness_width = roughness_width
        self._roughness = None
        self.build()

    @property
    def field_strength(self):
        return self._field_strength

    @field_strength.setter
    def field_strength(self, c):
        self._field_strength = c
        self.build()

    @property
    def scale(self):
        """Multipeak Mt. Fuji phenotypes without noise."""
        hd = np.empty([len(self.peaks), len(self.genotypes)])
        for i, peak in enumerate(self.peaks):
            newrow = self.hamming[i] * self.field_strength
            hd[i] = newrow

        min_hd = hd.min(0)  # Column-wise minimum value of array.
        max_min = np.amax(min_hd)  # Get the maximum value of the array for normalization.
        self._scale = 1 - (min_hd / max_min)  # Subtract from one -> Larger hamming dist. from peak = lower phenotype.
        return self._scale

    def build(self):
        """Construct phenotypes using a rough Mount Fuji model."""
        self.data.phenotypes = self.roughness + self.scale
