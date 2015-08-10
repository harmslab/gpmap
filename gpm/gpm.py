# Main mapping object to be used the epistasis models in this package.
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Outside imports
# ----------------------------------------------------------

import numpy as np

# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

# import different maps into this module
from gpm.base import BaseMap
from gpm.binary import BinaryMap
from gpm.mutations import MutationMap
from gpm.graph import Graph

# import utils used into module.
from gpm.utils import hamming_distance, binary_mutations_map, farthest_genotype, encode_mutations, construct_genotypes


class GenoPhenoMap(BaseMap):
    
    def __init__(self, wildtype, genotypes, phenotypes, errors=None, log_transform=False, mutations=None):
        """
            Construct a full genotype phenotype mapping object.
            
            Things included:
            ---------------
            1. Binary representation of all genotypes mapped to proper phenotypes
            2. Mutations mapped to their binary encoding
            3. NetworkX graph representation.
            
            Args:
            ----
            genotypes: array-like 
                list of all genotypes in system. Must be a complete system.
            phenotypes: array-like
                List of phenotypes in the same order as genotypes. 
            log_transform: boolean (default = False)
                Set to True to log tranform the phenotypes.
            mutations: dict
                Dictionary that maps each site indice to their possible substitution alphabet.
                
            Returns:
            -------
            GenoPhenoMap object
            
            
            GenoPhenoMap.Mutations --> includes all mutational mapping
            GenoPhenoMap.Binary --> genotypes mapped to their binary representations of the space. 
        
        """
        
        # Set mutations; if not given, assume binary space.
        if mutations is not None:
            self.mutations = mutations
        else:
            mutant = farthest_genotype(wildtype, genotypes)
            self.mutations = binary_mutations_map(wildtype, mutant)
        
        # Set initial properties fo GPM
        self.wildtype = wildtype
        self.genotypes = genotypes
        self.log_transform = log_transform
        self.phenotypes = phenotypes
        
        # Built the binary representation
        self._construct_binary()
        
        # If given errors, add them to map.
        if errors is not None:
            self.errors = errors
    
    @property
    def length(self):
        """ Get length of the genotypes. """
        return self._length    
    
    @property
    def n(self):
        """ Get number of genotypes, i.e. size of the system. """
        return self._n
    
    @property
    def log_transform(self):
        """ Boolean argument telling whether space is log transformed. """
        return self._log_transform
        
    @property
    def wildtype(self):
        """ Get reference genotypes for interactions. """
        return self._wildtype
        
    @property
    def mutations(self):
        """ Get the furthest genotype from the wildtype genotype. """
        return self._mutations

    @property
    def genotypes(self):
        """ Get the genotypes of the system. """
        return self._genotypes
        
    @property
    def phenotypes(self):
        """ Get the phenotypes of the system. """
        return self._phenotypes
    
    @property
    def errors(self):
        """ Get the phenotypes' errors in the system. """
        return self._errors

    @property    
    def indices(self):
        """ Return numpy array of genotypes position. """
        return self._indices
        
    # ----------------------------------------------------------
    # Getter methods for mapping objects
    # ----------------------------------------------------------   
    
    @property
    def geno2pheno(self):
        """ Return dict of genotypes mapped to phenotypes. """
        return self._map(self.genotypes, self.phenotypes)

    @property
    def geno2index(self):
        """ Return dict of genotypes mapped to their indices in transition matrix. """
        return self._map(self.genotypes, self.indices)
        
    @property
    def geno2binary(self):
        """ Return dictionary of genotypes mapped to their binary representation. """
        mapping = dict()
        for i in range(self.n):
            mapping[self.genotypes[self.Binary.indices[i]]] = self.Binary.genotypes[i] 
        return mapping
        
    # ----------------------------------------------------------
    # Setter methods
    # ----------------------------------------------------------
    
    @log_transform.setter
    def log_transform(self, boolean):
        """ True/False to log transform the space. """
        self._log_transform = boolean
    
    @genotypes.setter
    def genotypes(self, genotypes):
        """ Set genotypes from ordered list of sequences. """
        self._n = len(genotypes)
        self._length = len(genotypes[0])
        self._genotypes = np.array(genotypes)
        self._indices = np.arange(self.n)
        
    @wildtype.setter
    def wildtype(self, wildtype):
        """ Set the reference genotype among the mutants in the system. """
        self._wildtype = wildtype
        self.Mutations.wildtype = wildtype
    
    @mutations.setter
    def mutations(self, mutations):
        """ Set the mutation alphabet for all sites in wildtype genotype. 
         
            `mutations = { site_number : alphabet }`. If the site 
            alphabet is note included, the model will assume binary 
            between wildtype and derived.

            ``` 
            mutations = {
                0: [alphabet],
                1: [alphabet],

            }
            ```
        
        """
        if type(mutations) != dict:
            raise TypeError("mutations must be a dict")
        self._mutations = mutations
        self.Mutations = MutationMap()
        self.Mutations.mutations = mutations
        self.Mutations.n = len(mutations)

    @phenotypes.setter
    def phenotypes(self, phenotypes):
        """ Set phenotypes from ordered list of phenotypes 
            
            Args:
            -----
            phenotypes: array-like or dict
                if array-like, it musted be ordered by genotype; if dict,
                this method automatically orders the phenotypes into numpy
                array.
        """
        if type(phenotypes) is dict:
            self._phenotypes = self._if_dict(phenotypes)
        else:
            if len(phenotypes) != len(self._genotypes):
                raise ValueError("Number of phenotypes does not equal number of genotypes.")
            else:
                self._phenotypes = phenotypes

        # log transform if log_transform = True
        if self.log_transform is True:
            self._untransformed_phenotypes = self._phenotypes
            self._phenotypes = np.log10(self._phenotypes)

        
    @errors.setter
    def errors(self, errors):
        """ Set error from ordered list of phenotype error. 
            
            Args:
            -----
            error: array-like or dict
                if array-like, it musted be ordered by genotype; if dict,
                this method automatically orders the errors into numpy
                array.
        """
        # Order phenotype errors from geno2pheno_err dictionary
        if type(errors) is dict:
            errors = self._if_dict(errors)
        
        # For log-transformations of error, need to translate errors to center around 1,
        # then take the log.
        if self.log_transform is True:
            # Reference = http://onlinelibrary.wiley.com/doi/10.1002/sim.1525/epdf
            # \sigma_{log(f)}^{2} = log(1 + \sigma_{f}6{2}/mean(f)^{2}) 
            
            self._errors = np.array((   -np.sqrt(np.log10(1 + (errors**2)/self._untransformed_phenotypes**2)), 
                                        np.sqrt(np.log10(1 + (errors**2)/self._untransformed_phenotypes**2))))
                                                
            self.Binary._errors = np.array([self._errors[:,i] for i in self.Binary.indices]).T
        else:
            self._errors = errors
            self.Binary._errors = np.array([errors[i] for i in self.Binary.indices])
            
            
    # ------------------------------------------------------------
    # Displayed methods for mapping object
    # ------------------------------------------------------------
    
    def build_graph(self):
        """ Add a networkx DiGraph object representation of this space"""
        # Initialize network
        self.Graph = Graph(self)
         
    # ------------------------------------------------------------
    # Hidden methods for mapping object
    # ------------------------------------------------------------

    def _construct_binary(self):
        """ Encode the genotypes an ordered binary set of genotypes with 
            wildtype as reference state (ref is all zeros).
            
            This method maps each genotype to their binary representation
            relative to the 'wildtype' sequence.
        """
        self.Binary = BinaryMap()
        self.Binary.encoding = encode_mutations(self.wildtype, self.mutations)
        genotypes, self.Binary.genotypes = construct_genotypes(self.Binary.encoding)
        
        geno2index = self.geno2index
        self.Binary.indices = np.array([geno2index[genotypes[i]] for i in range(len(self.Binary.genotypes))])
        
        # Grab phenotypes if they exist. Otherwise, pass.
        try:
            self.Binary.phenotypes = np.array([geno2pheno[genotypes[i]] for i in range(len(self.Binary.genotypes))])
        except:
            pass
            
            