from functools import wraps
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from ..base import GenotypePhenotypeGraph
from .ellipses import draw_networkx_nodes_ellipses as draw_networkx_nodes
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
def labels(G, pos, ax, label_type="genotype", **kwargs):
    """Draw labels"""
    # Get labels from Graph.
    labels = dict([(node, G.node[node][label_type]) for node in G.nodes()])
    # Draw node labels
    nx.draw_networkx_labels(G, pos=pos, labels=labels, ax=ax,**kwargs)
    return ax

@checkG
def edges(G, pos, ax, **kwargs):
    """Draw edges.
    """
    nx.draw_networkx_edges(G, pos=pos, ax=ax, **kwargs)
    return ax

@checkG
def nodes(G, pos, ax, cmap='plasma', **kwargs):
    """Draw nodes
    """
    # Add color to nodes
    color = np.array([float(G.node[n]["phenotype"]) for n in G.nodes()])
    draw_networkx_nodes(G, pos, ax=ax, color=color, cmap=cmap, **kwargs)
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
def network(G, scale=1, vertical=True, figsize=(5,5), **kwargs):
    """Draw a generic plot for a genotype-phenotype map.
    """
    pos = positions.flattened(G, vertical=vertical, scale=scale)
    # separate keyword args.
    options = {"e":{}, "n":{}, "l":{}}
    for key, value in kwargs.items():
        options[key[0]][key[2:]] = value
    # init plot
    fig, ax = plt.subplots(figsize=figsize)
    ax = edges(G, pos, ax, **options["e"])
    ax = nodes(G, pos, ax, **options["n"])
    ax = labels(G, pos, ax, **options["l"])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.axis("off")
    return fig, ax, pos
