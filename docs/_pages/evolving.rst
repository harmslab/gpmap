Evolution
=========

The ``evolution`` module handles

Evolving Methods
----------------

Simulating evolution on genotype-phenotype maps can be done in a couple ways.
This page provides examples and tutorials that make this easy to learn. The
following functions can be imported fro teh ``gpmap.evolve`` package. All methods
must be given a source and target.

Monte Carlo
~~~~~~~~~~~

The ``monte_carlo`` function will use Monte Carlo sampling to walk through a
genotype-phentoype map. Given a source and target, it will return all steps taken
to move between the source and target.

.. code-block:: python

    from gpmap.evolve import monte_carlo


Evolutionary models
-------------------

GPMap comes with a few evolutionary models out of the box. They are live in the
``gpmap.evolve.models`` subpackage.

Gillespie Fixation
~~~~~~~~~~~~~~~~~~

The ``fixation`` function calculates a probability of fixing a new mutations given
the wildtype fitness and the new fitness using Gillespie's fixation model [1]_.

.. math::

    \pi_{i \rightarrow j} = \frac{1 - e^{-N \cdot \frac{f_j-f_i}{f_i}}}{1 - e^{-\frac{f_j-f_i}{f_i}}}

.. image:: ../_img/fixation.png
    :scale: 40 %
    :align: center


Source code:

.. code-block:: python

    Ns = [2, 5, 10, 100, ]
    fig, ax = plt.subplots(figsize=(5,3))

    fitness1 = 1
    fitness2 = np.linspace(-5,10,1000)
    sij = (fitness2 - fitness1)/abs(fitness1)

    for N in Ns:
        # Check the value of denominator
        denominator = 1 - np.exp(-N * sij)
        numerator = 1 - np.exp(- sij)
        # Calculate the fixation probability
        fixation = numerator / denominator
        ax.plot(sij, fixation, linewidth=3, label="N = " + str(N), alpha=.7)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xlabel("$\\frac{(f_1 -f_0)}{f_0}$", fontsize=16)
    ax.set_ylabel("$\pi$", fontsize=12)
    ax.set_title("Gillespie fixation probability")
    ax.legend(loc=2)



References
----------
.. [1] Gillespie, John H. Population genetics: a concise guide. JHU Press, 2010.
