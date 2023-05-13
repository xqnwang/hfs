## Grouped Variable Selection with Discrete Optimization
### Hussein Hazimeh, Rahul Mazumder, and Peter Radchenko

This is the accompanying code for our paper [Grouped Variable Selection with Discrete Optimization: Computational and Statistical Perspectives](https://arxiv.org/abs/2104.07084). 

This repo contains code for (i) approximate algorithms based on coordinate descent and combinatorial local search, and (ii) exact algorithms based on a custom branch-and-bound algorithm.

To get started please refer to [Demo.ipynb](https://github.com/hazimehh/L0Group)

## Prerequisites and Usage
The package is written in Python 3. It requires the following prerequisites:
- numpy
- scipy
- numba
- gurobi (only needed for the BnB algorithm)

See the Jupyter notebook Demo.ipynb for a demonstration on how to use the different algorithms.
