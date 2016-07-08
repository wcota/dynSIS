# Simulation of Markovian epidemic models on networks: SIS-II algorithm

## Synopsis

This code is a implementation of the SIS-II algorithm, as detailed in our paper (to be cited). It receives as input a network file, containing a list of edges and read, via terminal, the dynamical parameters.

## Dataset input

You need provide a file containing the list of edges (__in__ and __out__, two collumns). ID of the vertices must be enumerated sequentially as `1, 2, 3,..., N`, where `N` is the total number of vertices of the network. Here, we assume  __undirected__ and __unweighted__ networks without multiple neither self connections.

Consider, for example, a network with `N=5` vertices represented by:

```
1,2
1,3
2,4
2,5
3,4
```

Examples of datasets and their specifications will be available soon.

## Installation

In Linux and OSX, it is simple: just type ``make`` in the terminal in the *.f90 directory. If you need debugging, use ``make c=1``.

For Windows, however, you must compile all mod*.f90 files and the program code dynamics.f90. An example is:

```gfortran mod_read_tools.f90 mod_random.f90 mod_netdata.f90 dynamics.f90 -o dynamics```


## Use

Just run:

```./dynamics <edges_file> <output_file>```

``<output_file>`` will be written with the average of the fraction of infected vertices versus time.

## License

This code is under [GNU General Public License v3.0](http://choosealicense.com/licenses/gpl-3.0/).
