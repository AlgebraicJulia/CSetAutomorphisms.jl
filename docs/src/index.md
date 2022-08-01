# CSetAutomorphisms.jl

This library computes a canonical form for C-Sets, which are a data structure
akin to relational databases. The canonical form we want is one that is the same
no matter what order the database rows are in.

This problem has been extensively studied in the cases of graphs, which is a
C-set for a particular schema. However it is still very illustrative of how the
process works for C-sets in general. Note that, to efficiently represent the
graph on a computer, some order must be chosen for its vertices and edges. Yet,
we think of this as an implementation detail, not part of the 'real' graph. In
fact, we want to compare graphs (ask if they are equal or not) while ignoring
this order. We may wish to hash the graph to store it in a set or use it as a
key in a dictionary; we want this hash value to be independent of the order,
too. The computational problem is to search over all possible orderings for a
given graph and pick some *best* one to be our canonical representative.

This is described further in the following [blog post](https://www.algebraicjulia.org/blog/post/2022/01/cset-automorphisms/).

## Examples
[to do]