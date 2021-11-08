using BenchmarkTools
"""
This directory contains benchmarks for CSetAutomorphisms. To run all the
benchmarks, launch `julia --project=benchmark` and enter:

``` julia
using PkgBenchmark
import CSetAutomorphisms

benchmarkpkg(CSetAutomorphisms)
```

To run a specific set of benchmarks, use the `script` keyword argument, for
example:

``` julia
benchmarkpkg(CSetAutomorphisms; script="benchmark/FinSets.jl")
```
"""

const bench = SUITE = BenchmarkGroup()

using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams: TheoryWiringDiagram
using Random
include("../src/ColorRefine.jl")

# Color refine benchmarks
#########################

@acset_type WD(TheoryWiringDiagram)

# WDs
#----

# This could be generalized to a function making a 'random' C-set.
sizes = (Box=3, Wire=3, InPort=3, OutPort=3, OuterInPort=3, OuterOutPort=3,
         OutWire=3, InWire=3, PassWire=3)
wd = WD()
for (k,v) in pairs(sizes)
    add_parts!(wd, k, v)
end
for g in TheoryWiringDiagram.generators[:Hom]
    n, d_, cd_ = g.args
    d, cd = d_.args[1], cd_.args[1]
    set_subpart!(wd, n, rand(1:sizes[cd], sizes[d]))
end
[copy_parts!(wd, wd) for _ in 1:3] # add 8-fold symmetry
benchmark_crefine(bench, "color_refine_wd", color_refine(wd))

# Graphs
#-------
g1 = star_graph(Graph, 5)
g2 = path_graph(Graph, 5)
g = Graph(1)
[copy_parts!(g, h) for h in [g1, g2, g1, g2]]

benchmark_crefine(bench, "color_refine_graph", color_refine(g))

