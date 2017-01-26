Simulating
==========

The GPMap package comes with a suite of objects to simulate genotype-phenotype
maps following models in the literature. They are found in the ``gpmap.simulate``
module.

All Simulation objects inherit the ``GenotypePhenotypeMap`` object as their base
class. Thus, anything you can do with a GenotypePhenotypeMap, you can do with the
simulation objects.

NK model
--------

Construct a genotype-phenotype map using Kauffman's NK Model. [1]_
The NK fitness landscape is created using a table with binary, length-K,
sub-sequences mapped to random values. All genotypes are binary with length N.
The fitness of a genotype is constructed by summing the values of all
sub-sequences that make up the genotype using a sliding window across the full genotypes.

For example, imagine an NK simulation with :math:`N=5` and :math:`K=2`. To construct the fitness
for the 01011 genotype, select the following sub-sequences from an NK table:
"01", "10", "01", "11", "10". Sum their values together.

.. code-block:: python

    # import the NKSimulation class
    from gpmap.simulate import NKSimulation

    # Create an instance of the model. Using `from_length` makes this easy.
    gpm = NKSimulation.from_length(6, K=3)

House of Cards model
--------------------

Construct a 'House of Cards' fitness landscape. This is a limit of the NK model
where :math:`K=N`. It represents a fitness landscape with maximum roughness.


.. code-block:: python

    # import the HouseOfCardsSimulation class
    from gpmap.simulate import HouseOfCardsSimulation

    # Create an instance of the model. Using `from_length` makes this easy.
    gpm = HouseOfCardsSimulation.from_length(6)


Mount Fuji model
----------------

Construct a genotype-phenotype map from a Mount Fuji model. [2]_

A Mount Fuji sets a "global" fitness peak (max) on a single genotype in the space.
The fitness goes down as a function of hamming distance away from this genotype,
called a "fitness field". The strength (or scale) of this field is linear and
depends on the parameters `field_strength`.

Roughness can be added to the Mount Fuji model using a random `roughness` parameter.
This assigns a random roughness value to each genotype.

.. math::
    f(g) = \nu (g) - c \cdot d(g_0, g)

where :math:`\nu` is the roughness parameter, :math:`c` is the field strength, and :math:`d` is the
hamming distance between genotype :math:`g` and the reference genotype.

.. code-block:: python

    # import the HouseOfCardsSimulation class
    from gpmap.simulate import MountFujiSimulation

    # Create an instance of the model. Using `from_length` makes this easy.
    gpm = MountFujiSimulation.from_length(6)

    # add roughness, sampling from a range of values.
    gpm.set_roughness(range=(-1,1))

References
----------

.. [1] Kauffman, Stuart A., and Edward D. Weinberger. "The NK model of rugged fitness landscapes and its application to maturation of the immune response." Journal of theoretical biology 141.2 (1989): 211-245.
.. [2] Szendro, Ivan G., et al. "Quantitative analyses of empirical fitness landscapes." Journal of Statistical Mechanics: Theory and Experiment 2013.01 (2013): P01005.
