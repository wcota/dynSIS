# Simulation of Markovian epidemic models on networks: SIS-II algorithm

## Synopsis

What **is** it?

- [ ] to do

## To-do

- [ ] Write average value of samples
- [ ] Bash interface to read parameters
- [ ] Makefile
- [ ] Write informative text (reading network... running dynamics...)
- [ ] Write licence, URL, and DOI.

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

```ifort mod* dynamics.f90 -o dynamics```
or

```gfortran mod* dynamics.f90 -o dynamics```

- [ ] to do

## Use

```./dynamics <edges_file> <output_file>```

- [ ] to do

## License

This code is under [GNU General Public License v3.0](http://choosealicense.com/licenses/gpl-3.0/).
