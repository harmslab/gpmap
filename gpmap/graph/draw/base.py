from functools import wraps
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from ..base import GenotypePhenotypeGraph
from .ellipses import draw_networkx_nodes_ellipses
from . import positions

def checkG(func):
    """Check that all draw functions have GenotypePhenotypeGraph object as first
    argument.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        if type(args[0]) != GenotypePhenotypeGraph:
            raise TypeError("G must be an instance of GenotypePhenotypeGraph object.")
        return func(*args, **kwargs)
    return wrapper

@checkG
def labels(G, pos, ax, labels=None, label_type="genotype", **kwargs):
    """Draw labels"""
    # Get labels from Graph.
    if labels is None:
        labels = dict([(node, G.node[node][label_type]) for node in G.nodes()])
    # Draw node labels
    nx.draw_networkx_labels(G, pos=pos, labels=labels, ax=ax, **kwargs)
    return ax

@checkG
def edges(G, pos, ax, **kwargs):
    """Draw edges.
    """
    nx.draw_networkx_edges(G, pos=pos, ax=ax, **kwargs)
    return ax

@checkG
def nodes(G, pos, ax, cmap='plasma', color=None, colorbar=False, ellipses=True, **kwargs):
    """Draw nodes.
    """
    # Add color to nodes
    if color is None:
        color = np.array([float(G.node[n]["phenotype"]) for n in G.nodes()])
    if ellipses:
        ax, nodes = draw_networkx_nodes_ellipses(G, pos, ax=ax, color=color, cmap=cmap, **kwargs)
    else:
        nodes = nx.draw_networkx_nodes(G, pos, ax=ax, node_color=color, cmap=cmap, **kwargs)
    if colorbar is True:
        plt.colorbar(nodes, shrink=.3, aspect=5)
    return ax

@checkG
def arrow(G, pos, ax, source, target, scale=1, length=1,**kwargs):
    """Draw an arrow on graph
    """
    start = pos[source]
    end = pos[target]
    dx = length*(end[0] - start[0])
    dy = length*(end[1] - start[1])
    ax.quiver(start[0], start[1], dx, dy, scale_units="xy", angles="xy", scale=scale, **kwargs)
    return ax

@checkG
def path(G, pos, ax, path, scale=1, length=1, **kwargs):
    """Draw a path through network figure.
    """
    for i, node in enumerate(path[1:]):
        source = path[i]
        target = node
        arrow(G, pos, ax, source, target, scale=scale, length=length, **kwargs)
    return ax

@checkG
def network(G, scale=1, vertical=True, figsize=(5,5), ax=None, **kwargs):
    """Draw a standard networkx plot for a genotype-phenotype map. Flattens the
    network so that each row of nodes has the same number of mutations from
    wildtype.

    Parameters
    ----------
    G : GenotypePhenotypeGraph object
        genotype-phenotype network to plot
    scale : float
        scales the positions of the nodes in network

    Keyword Arguments
    -----------------
    Netork drawing is constructed from three separate networkx function calls:
    draw_networkx_edges, draw_networkx_nodes, and draw_networkx_node_labels. Keyword
    arguments are parsed and passed as kwargs to these separate functions. To direct
    these kwargs properly, edge kwargs should start with `e_`, node kwargs with
    `n_`, and label kwargs with `l_`.

    Returns
    -------
    fig : matplotlib Figure object
        Figure object.
    ax : matplotlib Axes object
        Subplot axes object
    pos : dictionary
        positions for each node.
    """
    pos = positions.flattened(G, vertical=vertical, scale=scale)
    # separate keyword args.
    options = {"e":{}, "n":{}, "l":{}}
    for key, value in kwargs.items():
        options[key[0]][key[2:]] = value
    # init plotx
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig=None
        #fig = ax.get_figure()

    ax = edges(G, pos, ax, **options["e"])
    ax = nodes(G, pos, ax, **options["n"])
    ax = labels(G, pos, ax, **options["l"])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.axis("off")
    return fig, ax, pos
