#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Outside imports
# ----------------------------------------------------------

import json
import pickle
import numpy as np
import pandas as pd

# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

# import different maps into this module
from . import mapping
from . import utils
from . import sample
from . import errors
from . import binary

DATA_KEYS = ("genotypes", "phenotypes", "stdeviations", "n_replicates")
OPTIONAL_KEYS = ("stdeviations", "n_replicates")

# ----------------------------------------------------------
# Exceptions
# ----------------------------------------------------------

class LoadingException(Exception):
    """Error when loading Genotype Phenotype map data. """

# ----------------------------------------------------------
# Base class for constructing a genotype-phenotype map
# ----------------------------------------------------------

class GenotypePhenotypeMap(mapping.BaseMap):
    """Main object for containing genotype-phenotype map data. Efficient memory storage,
    fast-ish mapping, graphing, and simulations.

    Parameters
    ----------
    wildtype : string
        wildtype sequence.
    genotypes : array-like
        list of all genotypes in system. Must be a complete system.
    phenotypes : array-like
        List of phenotypes in the same order as genotypes.
    mutations : dict
        Dictionary that maps each site indice to their possible substitution alphabet.
    n_replicates : int
        number of replicate measurements comprising the mean phenotypes
    include_binary : bool (default=True)
        Construct a binary representation of the space.
    """
    def __init__(self, wildtype, genotypes, phenotypes,
        stdeviations=None,
        mutations=None,
        n_replicates=1,
        include_binary=True,
        **kwargs):

        # Set mutations; if not given, assume binary space.
        if mutations is not None:
            # Make sure the keys in the mutations dict are integers, not strings.
            self.mutations = dict([(int(key), val) for key, val in mutations.items()])
        else:
            mutant = utils.farthest_genotype(wildtype, genotypes)
            mutations = utils.binary_mutations_map(wildtype, mutant)
            self.mutations = mutations

        # Set initial properties fo GPM
        self.wildtype = wildtype
        self.genotypes = genotypes
        self.phenotypes = phenotypes
        self.n_replicates = n_replicates
        self.stdeviations = stdeviations

        # Built the binary representation of the genotype-phenotype.
        # Constructs a complete sequence space and stores genotypes missing in the
        # data as an attribute, `missing_genotypes`.
        if include_binary:
            self.add_binary(self.wildtype)

        # Construct the error maps
        self._add_error()

    # ----------------------------------------------------------
    # Reading methods
    # ----------------------------------------------------------

    @classmethod
    def read_dataframe(cls, wildtype, dataframe, **kwargs):
        """Construct a GenotypePhenotypeMap from a dataframe."""
        # Required arguments
        df = DataFrames
        genotypes = df["genotypes"]
        phenotypes = df["phenotypes"]

        # Search for optional columns
        other_items = {}
        for key in OPTIONAL_KEYS:
            try: other_items[key] = df[key]
            except KeyError: pass

        # Initialize object.
        self = cls(wildtype, genotypes, phenotypes, **other_items)
        return self

    @classmethod
    def read_excel(cls, wildtype, fname, **kwargs):
        """"""
        df = pd.read_excel(fname)
        self = cls.read_dataframe(wildtype, df)
        return self

    @classmethod
    def read_csv(cls, fname, **kwargs):
        """"""
        df = pd.read_csv(fname)
        self = cls.read_dataframe(wildtype, df)
        return self

    @classmethod
    def read_json(cls, filename, **kwargs):
        """Load a genotype-phenotype map directly from a json file.
        The JSON metadata must include the following attributes

        Note
        ----
        Keyword arguments override input that is loaded from the JSON file.
        """
        # Open, json load, and close a json file
        f = open(filename, "r")
        data = json.load(f)
        f.close()

        # Grab all properties from data-structure
        necessary_args = ["wildtype", "genotypes", "phenotypes"]
        options = {
            "genotypes" : [],
            "phenotypes" : [],
            "wildtype" : [],
            "stdeviations": None,
            "mutations": None,
            "n_replicates": 1,
        }
        # Get all options for map and order them
        for key in options:
            # See if options are in json data
            try:
                options[key] = data[key]
            except:
                pass
        # Override any properties with manually entered kwargs passed directly into method
        options.update(kwargs)
        args = []
        for arg in necessary_args:
            val = options.pop(arg)
            args.append(val)
        # Create an instance
        gpm = cls(args[0], args[1], args[2], **options)
        return gpm

    # ----------------------------------------------------------
    # Writing methods
    # ----------------------------------------------------------

    def to_excel(self, filename, **kwargs):
        """Write genotype-phenotype map to excel spreadsheet.

        Keyword arguments are passed directly to Pandas dataframe to_excel method.

        Parameters
        ----------
        filename : str
            Name of file to write out.
        """
        self.df.to_excel(filename, **kwargs)

    def to_csv(self, filename):
        """Write genotype-phenotype map to csv spreadsheet.

        Keyword arguments are passed directly to Pandas dataframe to_csv method.

        Parameters
        ----------
        filename : str
            Name of file to write out.
        """
        self.df.to_csv(filename, **kwargs)

    def to_json(self, filename):
        """Write genotype-phenotype map to json file.
        """
        # Get metadata and make sure values are lists.
        data = self.metadata
        for key, val in data.items():
            if hasattr(val, "__iter__") and type(val) != dict:
                data[key] = list(val)
        # Write to file
        with open(filename, "w") as f:
            json.dump(data, f)

    # ----------------------------------------------------------
    # Properties methods
    # ----------------------------------------------------------

    @property
    def df(self):
        """Genotype-phenotype data in a DataFrame."""
        # Build dataframe
        data = {item : getattr(self, item) for item in DATA_KEYS}
        return pd.DataFrame(data, columns=DATA_KEYS)

    @property
    def length(self):
        """Get length of the genotypes. """
        return self._length

    @property
    def n(self):
        """Get number of genotypes, i.e. size of the genotype-phenotype map. """
        return self._n

    @property
    def wildtype(self):
        """Get reference genotypes for interactions. """
        return self._wildtype

    @property
    def mutant(self):
        """Get the farthest mutant in genotype-phenotype map."""
        _mutant = []
        _wt = self.wildtype
        for i in range(0,len(self.mutations)):
            site = _wt[i]
            options = self.mutations[i]
            if options is None:
                _mutant.append(_wt[i])
            else:
                for o in options:
                    if o != site:
                        _mutant.append(o)
        return "".join(_mutant)

    @property
    def mutations(self):
        """Get the furthest genotype from the wildtype genotype."""
        return self._mutations

    @property
    def genotypes(self):
        """Get the genotypes of the system."""
        return self._genotypes

    @property
    def missing_genotypes(self):
        """Genotypes that are missing from the complete genotype-to-phenotype map."""
        return self._missing_genotypes

    @property
    def complete_genotypes(self):
        """Array of sorted genotypes for the complete genotype space encoded by
        the mutations dictionary.

        **NOTE** Can only be set by the BinaryMap object.
        """
        try:
            return self._complete_genotypes
        except AttributeError:
            raise AttributeError("Looks like a BinaryMap has not been built yet for this map. Do this before asking for the complete_genotypes.")

    @property
    def phenotypes(self):
        """Get the phenotypes of the system. """
        return self._phenotypes

    @property
    def stdeviations(self):
        """Get stdeviations"""
        return self._stdeviations

    @property
    def n_replicates(self):
        """Return the number of replicate measurements made of the phenotype"""
        return self._n_replicates

    @property
    def index(self):
        """Return numpy array of genotypes position. """
        return self.genotypes.index

    # ----------------------------------------------------------
    # Setter methods
    # ----------------------------------------------------------

    @genotypes.setter
    def genotypes(self, genotypes):
        """Set genotypes from ordered list of sequences. """
        self._n = len(genotypes)
        self._length = len(genotypes[0])
        self._genotypes = pd.Series(genotypes)

    @wildtype.setter
    def wildtype(self, wildtype):
        """Set the reference genotype among the mutants in the system. """
        self._wildtype = wildtype

    @mutations.setter
    def mutations(self, mutations):
        """ Set the mutation alphabet for all sites in wildtype genotype.

        Examples
        --------
        `mutations = { site_number : alphabet }``. If the site
        alphabet is note included, the model will assume binary
        between wildtype and derived::

            mutations = {
                0: [alphabet],
                1: [alphabet],

            }
        """
        if type(mutations) != dict:
            raise TypeError("mutations must be a dict")
        # make sure keys are ints
        _mutations = {}
        for key, val in mutations.items():
            _mutations[int(key)] = val
        self._mutations = _mutations

    @phenotypes.setter
    def phenotypes(self, phenotypes):
        """Set phenotypes from ordered list of phenotypes."""
        self._phenotypes = pd.Series(phenotypes, index=self.index)

    @stdeviations.setter
    def stdeviations(self, stdeviations):
        """set stdeviations to array"""
        self._stdeviations = pd.Series(stdeviations, index=self.index)

    @n_replicates.setter
    def n_replicates(self, n_replicates):
        """Set the number of replicate measurements taken of phenotypes"""
        self._n_replicates = pd.Series(n_replicates, index=self.index)

    # ------------------------------------------------------------
    # Hidden methods for mapping object
    # ------------------------------------------------------------

    def _add_error(self):
        """Store error maps"""
        self.std = errors.StandardDeviationMap(self)
        self.err = errors.StandardErrorMap(self)

    def add_binary(self, wildtype):
        """Add a BinaryMap to the GenotypePhenotypeMap. The wildtype determines
        the encoding pattern. Wildtype sites are represented as 0's. 
        """
        self.binary = binary.BinaryMap(self, wildtype)

    # ------------------------------------------------------------
    # Hidden methods for mapping object
    # ------------------------------------------------------------

    def sample(self, n_samples=1, genotypes=None, fraction=1.0, derived=True):
        """Generate artificial data sampled from phenotype and percent error.

        Parameters
        ----------
        n_samples : int
            Number of samples to take from space

        fraction : float
            fraction of space to sample.

        Returns
        -------
        samples :  Sample object
            returns this object with all stats on experiment
        """
        if genotypes is None:
            # make sure fraction is float between 0 and 1
            if fraction < 0 or fraction > 1:
                raise Exception("fraction is invalid.")
            # fractional length of space.
            frac_length = int(fraction * self.n)
            # random genotypes and phenotypes to sample
            random_indices = np.sort(np.random.choice(range(self.n), size=frac_length, replace=False))
            # If sample must include derived, set the last random_indice to self.n-1
            if derived:
                random_indices[-1] = self.n-1
        else:
            # Mapping from genotypes to indices
            mapping = self.map("genotypes", "indices")
            # Construct an array of genotype indices to sample
            random_indices = [mapping[g] for g in genotypes]

        # initialize arrays
        phenotypes = np.empty((len(random_indices), n_samples), dtype=float)
        genotypes = np.empty((len(random_indices), n_samples), dtype='<U'+str(self.length))

        # If errors are present, sample from error distribution
        try:
            # Iterate through "seen" genotypes and sample from their distributions
            for i in range(len(random_indices)):
                index = random_indices[i]
                seq = self.genotypes[index]
                # Build genotype array
                genotypes[i] = np.array([seq for j in range(n_samples)])
                stdevs = self.err.upper
                phenotypes[i] = stdevs[index] * np.random.randn(n_samples) + self.phenotypes[index]
        except AttributeError:
            # Can't sample if no error distribution is given.
            if n_samples != 1:
                raise Exception("Won't create samples if sample error is not given.")
            genotypes = np.array([self.genotypes[i] for i in random_indices])
            phenotypes = np.array([self.phenotypes[i] for i in random_indices])

        # Create a sample object
        samples = sample.Sample(self, genotypes, phenotypes, random_indices)
        return samples

    def subspace(self, genotype1, genotype2=None, mutations=None):
        """Creates a sub-GenotypePhenotypeMap object between two genotypes.

        Parameters
        ----------
        genotype1 :
            wildtype genotype
        genotype2 :
            Farthest mutant sequence. Will assume from mutations dictionary if None
        mutations : dict
            Mutations dictionary

        Returns
        -------
        GenotypePhenotypeMap: obj
            Sub-GenotypePhenotypeMap between two sequences.
        """

        if genotype2 != None and len(genotype1) != len(genotype2):
            err = "genotypes must have the same length to calculate subspace.\n"
            raise ValueError(err)

        # Construct the mutations dictionary if not given (assumes binary)
        if mutations is None:
            mutations = utils.binary_mutations_map(genotype1, genotype2)
        # Construct binary encoding
        encoding = utils.encode_mutations(genotype1, mutations)
        # Construct the subspace
        wildtype = genotype1
        genotypes, binary = utils.construct_genotypes(encoding)
        # Get old genotype-phenotype mapping
        mapping = self.map("genotypes", "phenotypes")
        known, phenotypes = [], []
        for g in genotypes:
            try:
                phenotypes.append(mapping[g])
                known.append(g)
            except KeyError: pass
        # Get stdeviations
        if self.stdeviations is not None:
            mappingstd = self.map("genotypes", "stdeviations")
            stdeviations = [mappingstd[g] for g in known]
        else:
            stdeviations = None
        # Create GenotypePhenotypeMap object
        return GenotypePhenotypeMap(wildtype,
            known,
            phenotypes,
            stdeviations=stdeviations,
            log_transform=self.log_transform,
            mutations=mutations,
            n_replicates=self.n_replicates,
            logbase=self.logbase)
