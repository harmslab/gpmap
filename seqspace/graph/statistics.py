import numpy as np
import scipy as sp

def edge_flux(edges, transition_matrix):
    """Calculate the flux through set of edges from elements in transition
    matrix.

    Parameters
    ----------
    edges : list of tuples
        edges from a network graph.
    transition_matrix : 2d array
        Transition matrix. elements represent the probability of transitions from
        state i --> j and vice versa.

    Returns
    -------
    flux_matrix : 2d array
        Flux matrix calculated as state_frequencies * transition matrix
    """
    ### Calculate the state frequecies ###
    # Eigenvalues and Eigenvectors of transition matrix
    vals, vec = np.linalg.eig(transition_matrix)
    # Find the eigenvalue that == 1
    index = list(vals).index(1)
    state_freq = vec[:,index]
    ### Calculate the flux matrix ###
    flux_matrix = np.multiply(transition_matrix, state_freq)
    return flux_matrix / flux_matrix.sum(axis=1)


def current(edges, transition_matrix):
    """Calculate the flux through set of edges from elements in transition
    matrix.

    Parameters
    ----------
    edges : list of tuples
        edges from a network graph.
    transition_matrix : 2d array
        Transition matrix. elements represent the probability of transitions from
        state i --> j and vice versa.

    Returns
    -------
    flux_matrix : 2d array
        Flux matrix calculated as state_frequencies * transition matrix
    """
    ### Calculate the state frequecies ###
    # Eigenvalues and Eigenvectors of transition matrix
    vals, vl, vr = sp.linalg.eig(transition_matrix, left=True)
    # Find the eigenvalue that == 1
    index = list(vals).index(1)
    state_freq = vl[:,index]

    committor_plus = np.linalg.eig


    ### Calculate the flux matrix ###
    flux_matrix = np.multiply(transition_matrix, state_freq)
    return flux_matrix / flux_matrix.sum(axis=1)
