# GPM -- Python mapping utilities for constructing genotype-phenotype maps 

This module provides a mapping object/datastructures for genotype-phenotype maps. Easily converts these maps into Networkx objects and JSON datastructures. 

## Installation 

### Developers

Git must be installed to clone and contribute to this project.

1. Fork this repository on Github
2. Clone that repository locally
```
git clone <repo-url>
```
3. Navigate to this directory, and install (softly) this python package with 
```
cd <repo-name>
python setup.py develop
```
4. Add another remote link to the master version, call it `upstream`.
```
git remote add upstream <master-url-on-github>
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
git clone <master-url>
```

Navigate to this directory, and install this python package with 

```
python setup.py install
```

## Dependencies

The following modules are required for this to work. Also, the examples/tutorials are written in Jupyter notebooks and require IPython to be install. 

* [NetworkX](https://networkx.github.io/)
* [Numpy](http://www.numpy.org/)