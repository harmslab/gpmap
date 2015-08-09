# GPM -- Python mapping utilities for constructing genotype-phenotype maps 

This module provides mapping object/data-structures for managing genotype-phenotype data. 

1. Easily construct sequence spaces from lists of genotypes
2. Convert any set of sequences into binary representations for modeling.
3. Flexible framework for managing the type of mutations that occur in map.
4. Paired with LatticeGPM, can easily construct lattice protein sequence spaces.
5. Uses the power of [NetworkX](https://networkx.github.io) for network construction, enabling the use of their network algorithms.
6. Easily convert between Networkx objects and JSON data-structures.
7. Visualize networks with [NetworkViewer](https://github.com/harmslab/NetworkViewer) application. 

**NOTE**: Currently, these maps only work with complete spaces. We'll soon be working to make this more general.

## Installation 

### Developers

Git must be installed to clone and contribute to this project.

1. Fork this repository on Github
2. Clone that repository locally
```
git clone https://github.com/Zsailer/gpm
```
3. Navigate to this directory, and install (softly) this python package with 
```
cd gpm
python setup.py develop
```
4. Add another remote link to the master version, call it `upstream`.
```
git remote add upstream https://github.com/Zsailer/gpm
```
5. Start a branch locally from local master
```
git checkout -B <branch-name>
```
6. Make changes and commit to that branch.
```
git commit -a -m "<commit message>"
```
7. Push to your fork on github (which you called `upstream`).
```
git push upstream <branch-name>
```
8. Pull request the branch on Github into this master repository on Github.
## Users

Clone this repo locally:

```
git clone https://github.com/Zsailer/gpm
```

Navigate to this directory, and install this python package with 

```
python setup.py install
```

## Dependencies

The following modules are required for this to work. Also, the examples/tutorials are written in Jupyter notebooks and require IPython to be install. 

* [NetworkX](https://networkx.github.io/)
* [Numpy](http://www.numpy.org/)