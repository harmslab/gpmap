Network Graph
=============

``GenotypePhenotypeMap`` objects can be ported to NetworkX's Digraph object
using the ``add_networkx`` method. This appends DiGraph object to the ``Graph``
attribute of the GenotypePhenotypeMap. enables the use of all NetworkX's great
network analysis tools; see their docs_.

.. _docs: https://networkx.github.io/

Example
~~~~~~~

.. code-block:: python

    import networkx as nx

    gpm = GenotypePhenotypeMap.from_json("data.json")
    gpm.add_networkx()
    G = gpm.Graph

    nx.draw(G)

So, why didn't we start with NetworkX? Great question. Genotype-phenotype maps
scale horribly. NetworkX uses a dict-of-dict data structure, which is quite memory
intensive for larger genotype-phenotype maps. We store all values in numpy arrays.
If your data is small enough, you can easily port to a NetworkX Graph and take advantage
of all the functions provided by NetworkX.


Add evolutionary model
----------------------

Along with all the methods from

.. code-block:: python

    def adaptive(fitness1, fitness2):
        if fitness2 > fitness1:
            return 0
        else:
            return 1


    gpm.add_evolutionary_model(adaptive)


Interface
---------

.. autoclass:: gpmap.graph.base.GenotypePhenotypeGraph
    :members:
