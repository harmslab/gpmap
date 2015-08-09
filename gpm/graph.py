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
        for j in range(n_sub):
            # Get copy of reference string
            neighbor = str(reference[i])
            # Swap site for new string
            neighbor[i] = mutations[i][j]
            # Create a tuple of pair and append to list
            if mutation_label:
                # Add mutation
                mutation = reference[i] + str(i) + mutations[i][j]
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
        
        nodes = gpm.genotypes
        phenotypes = gpm.phenotypes
        reference = gpm.wildtype
        mutations = gpm.mutations
        errors = gpm.errors
        geno2binary = gpm.geno2binary
        
        for i in range(len(nodes)):
            # Add node to DiGraph
            self.add_node(
                nodes[i], 
                phenotype=phenotypes[i], 
                binary=geno2binary[genotypes[i]], 
                errors=errors[i]
            )
                   
            # Add edges from this node to Digraph
            edges = binary_neighbors(nodes[i], mutations, mutation_label=True)
            self.add_edges_from(edges)