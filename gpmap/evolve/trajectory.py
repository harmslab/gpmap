#
# class for holding trajectory information
#

from gpmap.utils import find_differences

class Trajectory(object):

    def __init__(self, gpm, trajectory, n=1):
        """Store useful information about a given trajectory.

        Parameters
        ----------
        gpm : GenotypePhenotypeMap object
            Genotype phenotype map with data for trajectory nodes
        trajectory : tuple
            Expects a tuple of genotypes in a trajectory
        n : int (default=1)
            number of times the trajectory was observed
        """
        self.gpm = gpm
        self.trajectory = trajectory
        self._counts = n

    @property
    def sites(self):
        """ Get a tuple with the sequence of sites mutated mutated. """
        return self.get_sites(self.trajectory)

    @property
    def mutations(self):
        """ Get a tuple of mutations in sequence. """
        return self.get_mutations(self.trajectory)

    @property
    def phenotypes(self):
        """ Get the phenotypes in a trajectory series. """
        return self.get_phenotypes(self.trajectory)

    @property
    def counts(self):
        """ Get the number of times a trajectory was visited. """
        return self._counts

    def add(self, n):
        """ add to count """
        self._counts += n

    def rm(self, n):
        """ rm from counts"""
        self._counts -= n

    @staticmethod
    def get_sites(trajectory):
        """ Get a tuple of the number of sites. """
        sites = tuple()

        # Iterate through trajectory and get the site number
        for i in range(1, len(trajectory)):
            index = find_differences(trajectory[i], trajectory[i-1])
            sites += (index[0]+1,)

        return sites

    @staticmethod
    def get_mutations(trajectory):
        """ Get a tuple of the number of sites. """
        mutations = tuple()

        # Iterate through trajectory and get the mutation
        for i in range(1, len(trajectory)):
            index = find_differences(trajectory[i], trajectory[i-1])
            mutations += (trajectory[i][index[0]],)

        return mutations

    def get_phenotypes(self, trajectory, log_transform=False):
        """ Get the phenotypes of each genotype in trajectory. """

        # Construct a mapping of genotypes to phenotypes and handle log tranforms
        if self.gpm.log_transform:

            # If user wants log transform, do it.
            if log_transform:
                mapping = self.gpm.map("genotypes", "phenotypes")

            # Else get the raw un-transformed
            else:
                mapping = self.gpm.map("genotypes", "phenotypes")

        else:
            mapping = self.gpm.map("genotypes", "phenotypes")

        # Iterate through trajectory and get the mutation
        phenotypes = tuple([mapping[t] for t in trajectory])

        return phenotypes
