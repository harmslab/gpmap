import numpy
import matplotlib as mpl
import networkx as nx

def draw_networkx_nodes_ellipses(G, pos,
    nodelist=None,
    height=20,
    width=30,
    angle=0,
    color='r',
    shape='o',
    alpha=1.0,
    cmap=None,
    vmin=None,
    vmax=None,
    ax=None,
    linewidths=None,
    label=None,
    **kwds):
    """Draw the nodes of the graph G.
    This draws only the nodes of the graph G.

    Parameters
    ----------
    G : graph
       A networkx graph
    pos : dictionary
       A dictionary with nodes as keys and positions as values.
       Positions should be sequences of length 2.
    ax : Matplotlib Axes object, optional
       Draw the graph in the specified Matplotlib axes.
    nodelist : list, optional
       Draw only specified nodes (default G.nodes())
    node_height : scalar or array
       Height of ellipse nodes (default=300).  If an array is specified it must be the
       same length as nodelist.
    node_width : scalar or array
       Width of ellipse nodes (default=300).  If an array is specified it must be the
       same length as nodelist.
    node_width : scalar or array
       Angle of major axis of ellipse nodes (default=300).  If an array is specified it must be the
       same length as nodelist.
    node_color : color string, or array of floats
       Node color. Can be a single color format string (default='r'),
       or a  sequence of colors with the same length as nodelist.
       If numeric values are specified they will be mapped to
       colors using the cmap and vmin,vmax parameters.  See
       matplotlib.scatter for more details.
    node_shape :  string
       The shape of the node.  Specification is as matplotlib.scatter
       marker, one of 'so^>v<dph8' (default='o').
    alpha : float
       The node transparency (default=1.0)
    cmap : Matplotlib colormap
       Colormap for mapping intensities of nodes (default=None)
    vmin,vmax : floats
       Minimum and maximum for node colormap scaling (default=None)
    linewidths : [None | scalar | sequence]
       Line width of symbol border (default =1.0)
    label : [None| string]
       Label for legend
    Returns
    -------
    matplotlib.collections.EllipseCollection
        `EllipseCollection` of the nodes.
    """

    import collections

    try:
        import numpy
        import matplotlib.pyplot as plt
        import matplotlib as mpl

    except ImportError:
        raise ImportError("Matplotlib required for draw()")
    except RuntimeError:
        print("Matplotlib unable to open display")
        raise

    if ax is None:
        ax = plt.gca()

    if nodelist is None:
        nodelist = list(G)

    if not nodelist or len(nodelist) == 0:  # empty nodelist, no drawing
        return None

    try:
        xy = numpy.asarray([pos[v] for v in nodelist])
    except KeyError as e:
        raise nx.NetworkXError('Node %s has no position.'%e)
    except ValueError:
        raise nx.NetworkXError('Bad value in node positions.')

    if cmap is not None:
        cmap = mpl.cm.get_cmap(cmap)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    node_collection = mpl.collections.EllipseCollection(widths=width,
        heights=height,
        angles=0,
        offsets=numpy.array(xy),
        cmap=cmap,
        norm=norm,
        transOffset=ax.transData,
        linewidths=linewidths)

    if type(color) == str:
        node_collection.set_color(color)
    else:
        if isinstance(alpha, collections.Iterable):
            color = apply_alpha(color, alpha, nodelist, cmap, vmin, vmax)
            alpha = None
        node_collection.set_array( color)
    node_collection.set_label(label)
    node_collection.set_zorder(2)
    ax.add_collection(node_collection)

    return ax, node_collection
