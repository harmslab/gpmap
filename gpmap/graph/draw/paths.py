import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from ..paths import flux
from .base import network
from .positions import flattened

def diff(G1, G2, source, target, node_scale=150, edge_scale=10, figsize=[2.5,2.5], ax=None):
    """
    """
    # Sanity check
    #if G1.gpm.genotypes.all() != G2.gpm.genotypes.all():
    #    raise Exception("Two GenotypePhenotypeGraphs are not comparable. "
    #    "Make sure the nodes are the same in both networks. "
    #    "Maybe try sorting the genotypes in G2 to match G1?")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig=None
        #fig = ax.get_figure()
    # Calculate fluxes
    flux1 = flux(G1, source, target)
    flux2 = flux(G2, source, target)
    edge_diff = flux2 - flux1
    node_flux1 = np.zeros(len(list(G1.nodes())))
    node_flux2 = np.zeros(len(list(G2.nodes())))
    for i, edge in enumerate(G1.edges()):
        # Get receiving node
        node = edge[1]
        if node not in [source, target]:
            node_flux1[node] += flux1[i]
            node_flux2[node] += flux2[i]
    node_diff = node_flux2 - node_flux1

    # Get node subsets
    node_increase = list(np.where(node_diff > 0)[0])
    node_decrease = list(np.where(node_diff < 0)[0])

    # Get edge subsets
    edge_increase = list(np.where(edge_diff > 0)[0])
    edge_decrease = list(np.where(edge_diff < 0)[0])
    edges = list(G1.edges())

    # Calculate positions
    positions = flattened(G1, vertical=True, scale=1)

    # Prepare plot
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.axis("off")

    # Plot increased nodes
    junk = nx.draw_networkx_nodes(G1, positions, ax=ax,
        nodelist=node_increase,
        node_size=node_scale * node_flux2[node_increase],
        node_color="b",
        linewidths=0,
    )
    # Overlay white center
    junk = nx.draw_networkx_nodes(G1, positions, ax=ax,
        nodelist=node_increase,
        node_size=node_scale * node_flux1[node_increase],
        node_color="w",
        linewidths=0,
    )
    # Plot decreased nodes
    junk = nx.draw_networkx_nodes(G1, positions, ax=ax,
        nodelist=node_decrease,
        node_size=node_scale * node_flux1[node_decrease],
        node_color="r",
        linewidths=0,
    )
    # Overlay white center
    junk = nx.draw_networkx_nodes(G1, positions, ax=ax,
        nodelist=node_decrease,
        node_size=node_scale * node_flux2[node_decrease],
        node_color="w",
        linewidths=0,
    )
    # Draw source and target
    junk = nx.draw_networkx_nodes(G1, positions, ax=ax,
        nodelist=[source,target],
        node_size=node_scale * .5,
        node_color="gray",
        linewidths=0,
    )
    # plot edge increase
    junk = nx.draw_networkx_edges(G1, positions, ax=ax,
        edgelist=[edges[i] for i in edge_increase],
        edge_color="b",
        width=edge_scale * abs(edge_diff[edge_increase]),
        arrows=False
    )
    # Plot edge decreases
    junk = nx.draw_networkx_edges(G1, positions, ax=ax,
        edgelist=[edges[i] for i in edge_decrease],
        edge_color="r",
        width=edge_scale * abs(edge_diff[edge_decrease]),
        arrows=False
    )
    return fig, ax


def edge_flux(G, source, target, width_scale=10, ax=None, **settings):
    """
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig=None
        #fig = ax.get_figure()
    # Calculate fluxes on edges
    fluxes = flux(G, source, target)

    # Settings to the network
    defaults = dict(
        figsize=[2.5,2.5],
        n_height=9,
        n_width=9,
        n_linewidths=0,
        n_vmax=round(max(G.gpm.phenotypes),2),
        n_vmin=round(max(G.gpm.phenotypes),2),
        n_colorbar=True,
        e_arrows=False,
        l_alpha=0.0,
    )
    defaults.update(**settings)
    fig, ax, pos = network(G, e_width = width_scale*fluxes, ax=ax, **settings)
    return fig, ax
