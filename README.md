# Optimized Gillespie algorithms for the simulation of Markovian epidemic processes on large and heterogeneous networks: SIS-OGA

This code is part of the article "[Optimized Gillespie algorithms for the simulation of Markovian epidemic processes on large and heterogeneous networks](https://doi.org/10.1016/j.cpc.2017.06.007)" [[ArXiv](https://arxiv.org/abs/1704.01557)].

[![license](https://img.shields.io/badge/licence-GPLv3-brightgreen.svg)](http://choosealicense.com/licenses/gpl-3.0/)
[![language](https://img.shields.io/badge/built%20with-Fortran-blue.svg)](https://gcc.gnu.org/fortran/)

### Fortran implementation

## Versions

(this) [Fortran implementation - for performance](https://github.com/wcota/dynSIS)

[Python implementation - learn and use](https://github.com/wcota/dynSIS-py)

[NetworkX Python implementation - range of options](https://github.com/wcota/dynSIS-networkx)

[GA Fortran implementation - Statistically exact, but NOT optimized](https://github.com/wcota/dynSIS-GA)

## Citation

Full bibliographic details: Computer Physics Communications 219C (2017) pp. 303-312

DOI information: 10.1016/j.cpc.2017.06.007

```
@article{COTA2017303,
title = "Optimized Gillespie algorithms for the simulation of Markovian epidemic processes on large and heterogeneous networks",
journal = "Computer Physics Communications",
volume = "219",
number = "",
pages = "303 - 312",
year = "2017",
note = "",
issn = "0010-4655",
doi = "http://dx.doi.org/10.1016/j.cpc.2017.06.007",
url = "http://www.sciencedirect.com/science/article/pii/S0010465517301893",
author = "Wesley Cota and Silvio C. Ferreira",
keywords = "Complex networks",
keywords = "Markovian epidemic processes",
keywords = "Gillespie algorithm",
abstract = "Numerical simulation of continuous-time Markovian processes is an essential and widely applied tool in the investigation of epidemic spreading on complex networks. Due to the high heterogeneity of the connectivity structure through which epidemic is transmitted, efficient and accurate implementations of generic epidemic processes are not trivial and deviations from statistically exact prescriptions can lead to uncontrolled biases. Based on the Gillespie algorithm (GA), in which only steps that change the state are considered, we develop numerical recipes and describe their computer implementations for statistically exact and computationally efficient simulations of generic Markovian epidemic processes aiming at highly heterogeneous and large networks. The central point of the recipes investigated here is to include phantom processes, that do not change the states but do count for time increments. We compare the efficiencies for the susceptible–infected–susceptible, contact process and susceptible–infected–recovered models, that are particular cases of a generic model considered here. We numerically confirm that the simulation outcomes of the optimized algorithms are statistically indistinguishable from the original GA and can be several orders of magnitude more efficient."
}
```

## Synopsis

This code is a implementation of the SIS-OGA (Optimized Gillespie Algorithm), as detailed in our [paper](https://doi.org/10.1016/j.cpc.2017.06.007). It receives as input a network file, containing a list of edges and read, via terminal, the dynamical parameters.

## Dataset input

You need to provide a file containing the list of edges (__in__ and __out__, two collumns). ID of the vertices must be enumerated sequentially as `1, 2, 3,..., N`, where `N` is the total number of vertices of the network. Here, we assume  __undirected__ and __unweighted__ networks without multiple neither self connections.

Consider, for example, a network with `N=5` vertices represented by:

```
1,2
1,3
2,4
2,5
3,4
```

Examples of datasets and their specifications are available at https://wcota.me/dynSISdatasets.

## Installation

In Linux and OSX, it is simple: just type ``make`` in the terminal in the *.f90 directory. If you need debugging, use ``make c=1``.

For Windows, however, you must compile all mod*.f90 files and the program code dynamics.f90. An example is:

```gfortran mod_read_tools.f90 mod_random.f90 mod_netdata.f90 dynamics.f90 -o dynamics```


## Use

If you want to manually input the dynamical parameters, just type:

```./dynamics <edges_file> <output_file>```

where ``<output_file>`` will be written with the average density of infected vertices versus time.

Alternatively, use (Linux):

```bash run.sh <edges_file> <output_file> <number of samples> <infection rate lambda> <maximum time steps> <fraction of infected vertices (initial condition)>```

_Example:_

```bash run.sh edges/s01.edges.dat "s01.lb0.002_100-samples.dat" 100 0.002 1000000 0.5```

## License

This code is under [GNU General Public License v3.0](http://choosealicense.com/licenses/gpl-3.0/).
