import numpy as np
import networkx as nx
from scipy.misc import comb
from matplotlib.cm import gray, spring
import matplotlib.pyplot as plt


class GraphDrawing(object):

    def __init__(self, graph, figsize=(6,6)):
        self._graph = graph
        self._configs = DrawingConfigs(self._graph, figsize=figsize)

    @property
    def configs(self):
        return self._config.__dict__

    def change(self, **kwargs):
        """ Change one of the configurations"""
        for kw in kwargs:
            setattr(self._config, kw, kwargs[kw])

        # Return fig, ax and config
        try:
            return self.fig, self.ax, config
        except:
            return self.configs

    def _draw_nodes(self):

        nx.draw_networkx_nodes(self._graph,
            pos=self._configs.pos,
            node_size=self._configs.node_size,
            node_color=self._configs.node_color,
            node_shape=self._configs.node_shape,
            alpha=self._configs.alpha,
            cmap=self._configs.cmap,
            vmin=self._configs.vmin,
            vmax=self._configs.vmax,

        )

    def flattened(self):
        pass

    def spring(self):
        pass




class DrawingConfigs(object):

    def __init__(self, graph, **kwargs):
        self._graph = graph
        self.pos = None
        self.alpha = 0.8
        self.cmap = gray
        self.node_color = colors
        self.node_size = 400
        self.arrows = False
        self.with_labels = True
        self.vmin = min(self._graph.gpm.phenotypes)
        self.vmax = max(self._graph.gpm.phenotypes)
        self.width = 0.5
        self.alpha = 1

        # Set any attributes desired.
        for kw in kwargs:
            setattr(self, kw, kwargs[kw])

    def spring_positions(self):
        """ """
        pos = nx.spring_layout(self._graph, iterations=300)
        return pos

    def flattened_positions(self, scale=1, vertical=False):
        """
            Get flattened positions for a genotype-phenotype graph.

            Input:
            -----
            space: GenotypePhenotypeGraph object

            Returns:
            -------
            positions: dict
                positions of all nodes in network (i.e. {index: [x,y]})
        """

        # Get the binary genotypes from GPM
        nodes = self.graph.gpm.binary.genotypes

        # Get mapping of binary genotypes to their graph indices
        mapping = self.graph.gpm.binary.map("genotypes", "indices")

        # Build an offset dictionary as we go...
        offsets = dict([(j, 0) for j in range(space.length+1)])

        # Init main positions dict
        positions = {}

        for i in range(space.n):
            # count the number of mutations for horizontal axis
            x = nodes[i].count("1")

            # Number of vertical positions
            pascal = comb(space.length, x)

            # Calculate the y position
            y = -(pascal-1) / 2 + offsets[x]

            # Add new position
            if vertical:
                positions[mapping[nodes[i]]] = [y, x*scale]
            else:
                positions[mapping[nodes[i]]] = [x*scale,y]

            # Iterate offset for that index on horizontal axis
            offsets[x] += 1.0

        self.pos = positions
        return positions

    def arrows(self):
        """ """
        for e in edges:

            arrows


        return



def draw_trajectories(G, trajectories, pos=None):
    pass

# ------------------------------------------------
# Methods for drawing trajectories on networks
# ------------------------------------------------

def edge_arrows(pos, edges):
    """ Maker a list of edge arrows. """
    arrows = list()
    for e in edges:
        arrows.append((pos[e[0]][0], pos[e[0]][1],
                       .9*(pos[e[1]][0]-pos[e[0]][0]),
                       .9*(pos[e[1]][1]- pos[e[0]][1]),
                       edges[e]))
    return arrows

def edge_weight(traj):
    """ Count the number of times each edge is visited. """
    edge_counter = dict()
    for t in traj:
        sequences = t
        edges = [(sequences[i-1],sequences[i]) for i in range(1,len(sequences))]
        for e in edges:
            if e in edge_counter:
                edge_counter[e] += traj[t]
            else:
                edge_counter[e] = traj[t]
    return edge_counter

def draw_space(G, pos=None):
    """ Draw the trajectories on Graph. """
    fig = plt.figure(figsize=[7,7])
    if pos is None:
        pos = nx.spring_layout(G, iterations=150)

    # Draw network
    colors = list()
    for n in G.nodes():
        colors.append(G.node[n]["phenotype"])

    nx.draw(G,pos, alpha=.8,
            cmap=gray,
            node_color=colors,
            node_size=400,
            arrows=False,
            with_labels=True,
            width=0.5,
            vmin = 0.94,
            vmax = 1.2,
           )
    return pos, fig

def draw_traj(G, traj, pos=None):
    """ Draw the trajectories on Graph. """
    fig = plt.figure(figsize=[7,7])
    if pos is None:
        pos = nx.spring_layout(G, iterations=150)

    # Draw network
    colors = list()
    for n in G.nodes():
        colors.append(G.node[n]["phenotype"])

    nx.draw(G,pos, alpha=.8,
            cmap=gray,
            node_color=colors,
            node_size=400,
            arrows=False,
            with_labels=True,
            width=2,
            vmin = 0.94,
            vmax = 1.2,
           )

    # Draw arrows
    edges = edge_weight(traj)
    arrows = edge_arrows(pos, edges)
    for a in arrows:
        plt.arrow(a[0], a[1], a[2], a[3], alpha=0.6, width=0.005*np.log(a[4]), head_width=0.05, head_length=0.05, fc='b', ec='k')
    return pos, fig
