# --------------------------------------------
# Module for running monte carlo simulations
# on genotype phenotype maps.
# --------------------------------------------

import numpy as np
from collections import Counter

class MaxIteractionError(Exception):
    """ Max number of iterations reached"""

class MonteCarloSimulation(object):
    
    
    def __init__(self, gpm, max_iter=1000):
        self.max_iter = max_iter
        self.gpm = gpm
        self._trajectories = Counter()
        self._n_trajectories = 0
        
    def run(self, n):
        """ Enumerate a given number of trajectories and add it to the trajectories property."""
        new_trajectories = self.enumerate_trajectories(self.gpm.Graph, n)
        self._trajectories += new_trajectories
        self._n_trajectories += n 
        
    @property
    def n_trajectories(self):
        """ Get the number of trajectories currently in the simulation. """    
        return self._n_trajectories
    
    @property
    def trajectories(self):
        """ Get the trajectories samples by simulation. """
        return self._trajectories

    @property
    def binary_trajectories(self):
        """ Get the trajectories with binary representation. """
        
        # Get the mapping between indices and binary genotypes
        index2binary = self.gpm.get_map("indices", "Binary.genotypes")
        
        mapping = {}
        # Iterate through trajectories and convert keys to binary repr.
        for key in self._trajectories:
            indices = list(key)
            sequences = tuple([index2binary[i] for i in indices])
            mapping[sequences] = self._trajectories[key]
            
        return mapping

    @staticmethod
    def trajectory(gpGraph, max_iter=1000):
        
        T = gpGraph.transition_matrix
        
        # Initial parameters for simulation
        iteration = 0            # Number of loop iterations
        finished = False         # Is the trajectory finished
        max_moves = len(T) - 1   # Make number of moves for forward trajectories only
        
        # Initialize the set of nodes visited in trajectory
        visited = (0,)
        
        # Begin simulation
        while finished is False and iteration < max_iter:
            
            # Separate transition probabilities into separate probabilty chucks to sample
            neighbor_probs = np.array(T[visited[-1]])[0]
            
            # Create a cumulative probability distribution from neighbor probs
            cumulative_prob_regions = np.array([sum(neighbor_probs[:i+1]) for i in range(len(neighbor_probs))])

            # Monte carlo time ... 
            mc_number = np.random.rand()    # Monte Carlo random number
            next_position = -1              # Tracker for next position
            new = -1                        # Trial moves (throwaway) for each monte carlo stel 
            
            # Monte Carlo sample the next move!
            while next_position == -1:
                new += 1
                if mc_number < cumulative_prob_regions[new]:
                    next_position = new
                    visited += (next_position,)
            
            if next_position == max_moves:
                finished = True
            
            # Add to iterations to avoid exceeding infinite looping simulations that get stuck!
            iteration += 1
            
        # Did we exceed the max iteration limit
        if iteration >= max_iter:
            raise MaxIterationsError("""Monte Carlo didn't finish. Reached max number of iterations!""")

        # Return the visited trajectory
        return visited
        
    @staticmethod
    def enumerate_trajectories(gpGraph, n):
        """ Enumerate n number of trajectories. 
        
            Returns a Counter object (type = dict) with key=trajectory, 
            and value=count of that trajectory. 
        """
        trajectories = [MonteCarloSimulation.trajectory(gpGraph) for i in range(n)]
        return Counter(trajectories)
        
    