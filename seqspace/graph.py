import numpy as np
from networkx import DiGraph

# --------------------------------------------------------
# Utils for network building
# --------------------------------------------------------

def binary_neighbors(reference, mutations, mutation_label=False):
    """ Return neighbors to reference string using mutations dictionary
        and return neighbor pairs.
    """
    neighbor_pairs = list()
    n_sites = len(reference)
    
    for i in range(n_sites):
        n_sub = len(mutations[i])
        possible = list(mutations[i])
        possible.remove(reference[i])
        # Create a tuple of pair and append to list
        for j in range(n_sub-1):
            neighbor = list(reference)
            neighbor[i] = possible[j]
            neighbor = "".join(neighbor)
            if mutation_label:
                # Add mutation
                mutation = reference[i] + str(i) + possible[j]
                neighbor_pairs.append((reference, neighbor, {"mutation":mutation}))
            else:    
                neighbor_pairs.append((reference, neighbor))
    return neighbor_pairs

# --------------------------------------------------------
# Utils for network building
# --------------------------------------------------------

class Graph(DiGraph):
    
    def __init__(self, gpm):
        """ Construct a DiGraph network from gpm. """
        
        # initialize the DiGraph object
        super(Graph, self).__init__()
        
        # Grab properties of parentmapping object
        nodes = gpm.genotypes
        phenotypes = gpm.phenotypes
        reference = gpm.wildtype
        mutations = gpm.mutations
        geno2binary = gpm.geno2binary
        
        for i in range(len(nodes)):
            # If no error is present, store None
            try:
                error = gpm.errors[i]
            except AttributeError:
                error = None
            
            # Add node to DiGraph
            self.add_node(
                nodes[i], 
                phenotype=phenotypes[i], 
                binary=geno2binary[nodes[i]], 
                errors=error
            )
                   
            # Add edges from this node to Digraph
            edges = binary_neighbors(nodes[i], mutations, mutation_label=True)
            self.add_edges_from(edges)