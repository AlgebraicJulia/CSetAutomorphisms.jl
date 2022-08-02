# CSetAutomorphisms.jl [![Documentation](https://github.com/kris-brown/ModelEnumeration.jl/workflows/Documentation/badge.svg)](https://AlgebraicJulia.github.io/CSetAutomorphisms.jl/dev/)
![Tests](https://github.com/AlgebraicJulia/CSetAutomorphisms.jl/workflows/Tests/badge.svg)

[Attributed C-sets](https://arxiv.org/pdf/2106.04703.pdf) encompass a broad class of data structures, including many generalizations of graphs (e.g. [directed](https://www.algebraicjulia.org/blog/post/2020/09/cset-graphs-1/), [symmetric](https://www.algebraicjulia.org/blog/post/2020/09/cset-graphs-2), [reflexive](https://www.algebraicjulia.org/blog/post/2021/04/cset-graphs-3/)), tabular data (e.g. [data frames](https://pandas.pydata.org/pandas-docs/stable/user_guide/dsintro.html)), and combinations of the two (e.g. weighted graphs, [relational databases](https://en.wikiversity.org/wiki/Relational_Databases/Introduction)). This repo generalizes the [Nauty](https://pallini.di.uniroma1.it/Introduction.html) algorithm, which produces canonical members of an isomorphism class, to cover any data structure encompassed by C-Sets.

More background is found in the documentation (click above) or this [blog](https://www.algebraicjulia.org/blog/post/2022/01/cset-automorphisms/) post.

## To do
- An upcoming refactor of Catlab using [CompTime.jl](https://github.com/olynch/CompTime.jl) will make it feasible to reimplement the core algorithm using code custom generated for a given C-Set, offering new opportunities for performance improvements.
- Rather than just computing a canonical hash of a C-set's isomorphism class, we want to compute and represent the full automorphism groups, which can be applied to solving certain problems (e.g. searching for homomorphisms from A to B, up to isomorphism).

