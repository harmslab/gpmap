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
        genotypes = gpm.genotypes
        phenotypes = gpm.phenotypes
        reference = gpm.wildtype
        mutations = gpm.mutations
        geno2binary = gpm.get_map("genotypes", "Binary.genotypes")
        geno2index = gpm.get_map("genotypes", "indices")
        
        for i in range(len(genotypes)):
            # If no error is present, store None
            try:
                error = float(gpm.errors[i])
            except AttributeError:
                error = None
            
            # Add node to DiGraph
            self.add_node(
                int(geno2index[genotypes[i]]),           # genotype index
                genotype=str(genotypes[i]),              # genotype
                binary=str(geno2binary[genotypes[i]]),   # binary representation
                phenotype=float(phenotypes[i]),            # phenotype
                value=float(phenotypes[i]),                    # same as phenotype
                errors=error                        # error in phenotype
            )
                   
            # Add edges from this node to Digraph
            mutation_label = True
            edges = binary_neighbors(genotypes[i], mutations, mutation_label=mutation_label)
            
            if mutation_label is True:
                edges_ = [(int(geno2index[e[0]]), int(geno2index[e[1]]), e[2]) for e in edges]

            self.add_edges_from(edges_)