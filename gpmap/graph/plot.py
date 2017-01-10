import matplotlib.pyplot as plt

def path_spectrum(probs, ax=None, figsize=(4,1), ybounds=None, **kwargs):
    """Plot a bar graph of all forward paths in G versus their relative probabilities.
    """
    # Plot on preconfigured axis
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    # Set ybounds
    if ybounds is None:
        ymax = max(probs)
        ymax =
        ybounds = [0, max()]
    # Plot
    ax.bar(range(len(probs)), probs, width=1, **kwargs)
    # Set spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_bounds(0, len(paths))
    ax.spines["left"].set_bounds(0,0.15)
    # Set ticks
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')




    ax.set_yticks(np.linspace(0,0.09, 4))
    ax.set_xticks(np.linspace(0,120, 5))
    # Axis size
    ax.axis([0, len(paths)+.5] + yaxis)
    # Draw a grid
    ax.grid(linestyle=":", color="k")
    ax.yaxis.grid(False)
    return fig, ax
