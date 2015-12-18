import numpy as np
import networkx as nx

# --------------------------------------------------------
# Utils for network building
# --------------------------------------------------------

def binary_neighbors(reference, mutations, mutation_labels=False):
    """ Return neighbors to reference string using mutations dictionary
        and return neighbor pairs.
        
        Returns:
        -------
        neighbor_pairs = [(genotype1, genotype2), ...]
                       or
                       = [(genotype1, genotype2, {mutations: "0A1"}), ...]
    """
    neighbor_pairs = list()
    n_sites = len(reference)

    for i in range(n_sites):
        
        # Calculate the number of mutations possible at a given sites
        if mutations[i] == None:
            n_sub = 1
        
        else:
            n_sub = len(mutations[i])
            
            # Remove the reference mutation
            possible = list(mutations[i])
            possible.remove(reference[i])
        
        # Create a tuple of pair and append to list
        for j in range(n_sub-1):
            # Switch out site i with each mutation possible at that site
            neighbor = list(reference)
            neighbor[i] = possible[j]
            neighbor = "".join(neighbor)
            
            # If user wants the mutation to be defined with a string (i.e. `A42G`), do that.
            if mutation_labels:
                # Add mutation
                mutation = reference[i] + str(i) + possible[j]
                neighbor_pairs.append((reference, neighbor, {"mutation":mutation}))
            
            else:
                neighbor_pairs.append((reference, neighbor, {}))
                
    return neighbor_pairs

# --------------------------------------------------------
# Utils for network building
# --------------------------------------------------------

class GenotypePhenotypeGraph(nx.DiGraph):

    def __init__(self, gpm):
        """ Construct a DiGraph network from gpm. """

        # initialize the DiGraph object
        super(GenotypePhenotypeGraph, self).__init__()

        self.gpm = gpm
            
    
    def add_gpm_node(self, index, genotype=None, binary=None, phenotype=None, value=None, errors=None, **kwargs):
        """ ADD node to networkx graph. """
        self.add_node(index, 
            genotype=genotype,
            binary=binary,
            phenotype=phenotype,
            value=value,
            errors=errors,
            **kwargs
        )

        
    def add_gpm_edges(self, ebunch, transition_func=None):
        """ Method for adding edges to the graph. """
        
        # Check whether the edges are genotypes or indices --> convert to indices.

        if type(ebunch[0][0]) is str or type(ebunch[0][0]) is np.str_:
            geno2index = self.gpm.get_map("genotypes", "indices")
            node = lambda x: geno2index[x]
        else:
            node = lambda x: x
    
        # Add edge with data (include transition data if given)
        edges = list
        for edge in ebunch:
            
            # Get indices of neighbor nodes
            index = node(edge[0])
            index2 = node(edge[1])
            attributes = edge[2]
            
            # Run transition calculation
            if transition_func is not None:            
                # Calculate the transition function 
                attributes["transition"] = transition_func(self.gpm.phenotypes[index], self.gpm.phenotypes[index2])
            
            self.add_edge(index, index2, attributes)
        
        
    def build_graph(self, transition_func=None, mutation_labels=False):
        """ Build the graph. """
        try:
            errors = np.array(self.gpm.err.upper)
        except:
            errors = [None for i in range(self.gpm.n)]

        for i in range(self.gpm.n):
            
            # If no error is present, store None
            geno2index = self.gpm.get_map("genotypes", "indices")

            # Construct nodes to add to the graph
            self.add_gpm_node(
                int(geno2index[self.gpm.genotypes[i]]),                     # genotype index
                genotype=str(self.gpm.genotypes[i]),            # genotype
                binary=str(self.gpm.Binary.genotypes[i]),      # binary representation
                phenotype=float(self.gpm.phenotypes[i]),        # phenotype
                value=float(self.gpm.phenotypes[i]),            # same as phenotype
                errors=errors[i]                                # error in phenotype
            )
        
            # Construct a set of edge labels to add to Graph
            edges = binary_neighbors(self.gpm.genotypes[i],
                        self.gpm.mutations,
                        mutation_labels=mutation_labels
            )

            # Add edges to map
            self.add_gpm_edges(edges, transition_func=transition_func)
            
    @property
    def transition_matrix(self):
        """ Get transition matrix of the graph. Only works if transitions is function is set."""
        return np.nan_to_num(
            nx.attr_matrix(
                self,
                edge_attr="transition",         # Get the transitions value
                normalized=True,                # Normalize the rows
                rc_order=self.gpm.indices       # Order the matrix
                )
        )