import numpy as np
import networkx as nx

# --------------------------------------------------------
# Utils for network building
# --------------------------------------------------------

def binary_neighbors(reference, mutations, mutation_labels=False):
    """ Return neighbors to reference string using mutations dictionary
    and return neighbor pairs.

    Returns
    -------
    neighbor_pairs : list of tuple pairs
        [(genotype1, genotype2), ...] or [(genotype1, genotype2, {mutations: "0A1"}), ...]
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
    """Construct a DiGraph network from gpm.
    """
    def __init__(self, gpm):
        # initialize the DiGraph object
        super(GenotypePhenotypeGraph, self).__init__()
        self.gpm = gpm
        self.built = False
        self.transition_model = None

    def add_gpm_node(self, index, genotype=None, binary=None, phenotype=None, value=None, errors=None, **kwargs):
        """Add node to networkx graph. """
        self.add_node(index,
            genotype=genotype,
            binary=binary,
            phenotype=phenotype,
            value=value,
            errors=errors,
            **kwargs
        )

    def add_gpm_edges(self, ebunch):
        """Method for adding edges to the graph. """
        # Check whether the edges are genotypes or indices --> convert to indices.
        if type(ebunch[0][0]) is str or type(ebunch[0][0]) is np.str_:
            genotypes = self.gpm.complete_genotypes
            geno2index = dict(zip(self.gpm.complete_genotypes, range(len(genotypes))))
            node = lambda x: geno2index[x]
        else:
            node = lambda x: x

        # Add edge with data (include transition data if given)
        edges = list
        for edge in ebunch:
            # Get indices of neighbor nodes
            try:
                index = node(edge[0])
                index2 = node(edge[1])
                attributes = edge[2]
                self.add_edge(index, index2, **attributes)
            except KeyError:
                pass

    def add_evolutionary_model(self, model, *args, **kwargs):
        """Add an evolutionary model to the genotype phenotype graph. The model
        argument describes the probability of fixation model for each edge in the graph.

        The main assumption of this method is that each node has an attribute named
        values, which is the fitness of that genotype.
        """
        self.transition_model = model
        for e in self.edges():
            # Acquire states.
            i = e[0]
            j = e[1]
            state_i = self.node[i]
            state_j = self.node[j]
            # Calculate the fixation probability for this age.
            fixation = model(state_i["value"], state_j["value"], *args, **kwargs)
            # Set the edge value
            self.edge[i][j]["fixation"] = fixation


    def _build(self, log_transform=False, include_missing=True, mutation_labels=False):
        """Attach a Network DiGraph to GenotypePhenotypeMap object."""
        # If log transform
        if log_transform:
            phenotypes = self.gpm.log.phenotypes
        else:
            phenotypes = self.gpm.phenotypes
        if self.gpm.stdeviations is None:
            errors = [None for i in range(self.gpm.n)]
        else:
            errors = np.array(self.gpm.err.upper)
        edges = []
        for i in range(self.gpm.n):
            # If no error is present, store None
            geno2index = self.gpm.map("genotypes", "indices")
            # Construct nodes to add to the graph
            self.add_gpm_node(
                int(geno2index[self.gpm.genotypes[i]]),         # genotype index
                genotype=str(self.gpm.genotypes[i]),            # genotype
                binary=str(self.gpm.binary.genotypes[i]),       # binary representation
                phenotype=phenotypes[i],                            # phenotype
                value=phenotypes[i],                                # same as phenotype
                errors=errors[i]                                # error in phenotype
            )
            # Construct a set of edge labels to add to Graph
            edges += binary_neighbors(self.gpm.genotypes[i],
                        self.gpm.mutations,
                        mutation_labels=mutation_labels)
        # Add missing nodes if true.
        if include_missing:
            genotypes = self.gpm.missing_genotypes
            for i, genotype in enumerate(genotypes):
                # Construct nodes to add to the graph
                index = int(self.gpm.n + i)
                self.add_gpm_node(
                    index, # genotype index
                    genotype=str(genotype),  # genotype
                    binary=str(self.gpm.binary.complete_genotypes[index]))
                # Construct a set of edge labels to add to Graph
                edges += binary_neighbors(genotype,
                            self.gpm.mutations,
                            mutation_labels=mutation_labels)
        # Add edges to map
        self.add_gpm_edges(edges)

    @property
    def transition_matrix(self):
        """Get transition matrix of the graph. Only works if transitions is
        function is set.

        Calculates each element as

        .. math:

            p_{i \rightarrow j} = 1 / N \cdot \pi_{i \rightarrow j}

        where N is the number of neighbors and pi is a fixation probability
        """
        if self.transition_model is None:
            raise Exception("""A transition model must be set. checkout `add_evolutionary_model` method.""")
        # Get a matrix of fixation probabilities.
        matrix = np.nan_to_num( nx.attr_matrix(self, edge_attr="fixation",rc_order=self.gpm.indices))
        # scale fixation probabilities by the condition probability.
        matrix = matrix / self.gpm.binary.length #+ 1
        # Calculate the self probabilities (diagonal)
        for i in range(len(matrix)):
            matrix[i,i] = 1 - matrix[i].sum()
        return matrix
