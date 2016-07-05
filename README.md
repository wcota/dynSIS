# Simulation of Markovian epidemic models on large networks: SIS algorithm

## Synopsis

What **is** it?

- [ ] to do

## Dataset input

All you need is a file containing the list of **unique** edges (__in__ and __out__, two collumns) between all nodes of the network. ID of each vertex must be larger or equal to `1`, and be sequential: `1, 2, 3,..., N`, where `N` is the total number of vertices of the network. Here, we assume only __undirected__ networks and edges weight equal to 1.

As an example, consider a network of `N=5` vertices. `1` is connected to `2,3`, `2` to `1,4,5`, `3` to `1,4`, `4` to `2,3` and `5` to `2`. So, the file would be:

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

- [ ] to do

## Use

```./dynamics <edges_file> <output_file>```

- [ ] to do

## License

This code is under [GNU General Public License v3.0](http://choosealicense.com/licenses/gpl-3.0/).
