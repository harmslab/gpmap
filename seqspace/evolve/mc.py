# --------------------------------------------
# Module for running monte carlo simulations
# on genotype phenotype maps.
# --------------------------------------------

import numpy as np
import networkx as nx
from seqspace.plotting import TrajectoriesPlotting
from collections import Counter, OrderedDict

class MaxIterationsError(Exception):
    """ Max number of iterations reached"""


class Trajectories(object):
    
    def __init__(self, gpm, source, target, trajectory_data={}):
        """ Object for holding trajectories from MC simulation. """ 
        self._gpm = gpm  
        self.source = source
        self.target = target 
        self.trajectory_data = trajectory_data
        self.Plot = TrajectoriesPlotting(self)
        
    @property
    def nodes(self):
        """ Get the trajectory data with node index as key"""
        return self.sort_dict(self.trajectory_data)

    @property
    def n(self):
        """ Get the total count of all trajectories. """
        return sum(list(self.nodes.values()))
        
    def _forward(self, trajectories):
        """ Return only the trajectories that move foward."""
        forward = [(key, trajectories[key]) for key in trajectories if len(key)-1 == self._gpm.Mutations.n]
        return OrderedDict(forward)
        
    @property
    def possible(self):
        """List the set of possible forward trajectories. """
        return [tuple(path) for path in nx.all_shortest_paths(self._gpm.Graph, source=self.source, target=self.target)]
         
    @property
    def forward(self):
        """ Get forward trajectories. """
        return self._forward(self.nodes)
        
    def _binary(self, trajectories):
        """ Replace the keys in the trajectories dictionary for binary genotypes. """
        # Get the mapping between indices and binary genotypes
        index2binary = self._gpm.get_map("indices", "Binary.genotypes")
        
        # New dictionary
        mapping = OrderedDict()
        
        # Iterate through trajectories and convert keys to binary repr.
        for key in trajectories:
            indices = list(key)
            sequences = tuple([index2binary[i] for i in indices])
            mapping[sequences] = trajectories[key]
        
        _mapping = self.sort_dict(mapping)
        return _mapping
        
    @property
    def binary(self):
        """ Get trajectory dict with binary genotypes as keys. """
        return self._binary(self.nodes)
        
    @property
    def fbinary(self):
        """ Get forward binary trajectories."""
        return self._binary(self.forward)    
        
    @staticmethod
    def sort_dict(dictionary):
        """ Sort a dictionary by its value"""
        sorted_counts = sorted(list(dictionary.values()), reverse=True)
        ordered = OrderedDict()
        for s in sorted_counts:
            for key in dictionary:
                if dictionary[key] == s:
                    ordered[key] = s
        return ordered

class MonteCarloSimulation(object):
    
    
    def __init__(self, gpm, source, target, max_iter=1000):
        self.max_iter = max_iter
        self.gpm = gpm
        self.source = source
        self.target = target
        self._trajectories = Counter()
        self.Trajectories = Trajectories(self.gpm, self.source, self.target, self._trajectories)
        
        
    def run(self, n):
        """ Enumerate a given number of trajectories and add it to the trajectories property."""
        new_trajectories = self.enumerate_trajectories(self.gpm.Graph, n, self.source, self.target, max_iter=self.max_iter)
        self._trajectories += new_trajectories

    @staticmethod
    def trajectory(gpGraph, source, target, max_iter=1000):
        
        T = gpGraph.transition_matrix
        
        # Initial parameters for simulation
        iteration = 0            # Number of loop iterations
        finished = False         # Is the trajectory finished
        #max_moves = len(T) - 1   # Make number of moves for forward trajectories only
        
        # Initialize the set of nodes visited in trajectory
        visited = (source,)
        
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
            while next_position == -1 and new + 1 < len(cumulative_prob_regions):
                new += 1

                if mc_number < cumulative_prob_regions[new]:
                    next_position = new
                    visited += (next_position,)
            
            #if next_position == max_moves:
            if next_position == target:
                
                finished = True
            
            # Add to iterations to avoid exceeding infinite looping simulations that get stuck!
            iteration += 1
            
        # Did we exceed the max iteration limit
        if iteration >= max_iter:
            print(visited)
            raise MaxIterationsError("Monte Carlo didn't finish. Reached max number (" + str(max_iter) + ") of iterations")

        # Return the visited trajectory
        return visited
        
    @staticmethod
    def enumerate_trajectories(gpGraph, n, source, target, max_iter=1000):
        """ Enumerate n number of trajectories. 
        
            Returns a Counter object (type = dict) with key=trajectory, 
            and value=count of that trajectory. 
        """
        trajectories = [MonteCarloSimulation.trajectory(gpGraph, source, target, max_iter) for i in range(n)]
        return Counter(trajectories)
        
    