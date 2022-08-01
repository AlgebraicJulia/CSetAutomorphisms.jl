using BenchmarkTools
const bench = SUITE = BenchmarkGroup()

using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams: TheoryWiringDiagram
using Catlab.Graphs
using Catlab.Theories
using Random: rand
using CSetAutomorphisms

# Color refine benchmarks
#########################
# @acset_type WD(TheoryWiringDiagram)
# WDs
#----
# This could be generalized to a function making a 'random' C-set.
# sizes = (Box=3, Wire=3, InPort=3, OutPort=3, OuterInPort=3, OuterOutPort=3,
#          OutWire=3, InWire=3, PassWire=3)
# wd = WD()
# for (k,v) in pairs(sizes)
#   add_parts!(wd, k, v)
# end
# for g in TheoryWiringDiagram.generators[:Hom]
#   n, d_, cd_ = g.args
#   d, cd = d_.args[1], cd_.args[1]
#   set_subpart!(wd, n, rand(1:sizes[cd], sizes[d]))
# end
# [copy_parts!(wd, wd) for _ in 1:3] # add 8-fold symmetry

# bench["Color saturate a random wiring diagram:"] = @benchmarkable color_saturate(wd)

# Graphs
#-------


# G = Graph(1)
# [copy_parts!(G, x) for x in [g1,g2]]
# canonical_hash(G); call_nauty(G)
# bench["Canonical form, nauty:"]= @benchmarkable call_nauty(G)
# bench["Canonical form, csets:"]= @benchmarkable canonical_hash(G)
n_reps = 10
times = []
for i in [3,5,7,9]
  println(i)
  push!(times,[])
  G = star_graph(ReflexiveGraph, i) ⊕ path_graph(ReflexiveGraph, 3)
  add_edge!(G, i, i+1)
  X = G ⊗ G
  X′ = elements(X)
  for _ in 1:n_reps
    e1 = @elapsed(canonical_hash(X))
    e2 = @elapsed(canonical_hash(X′))
    e3 = @elapsed(call_nauty(X))
    push!(times[end], (e1,e2,e3))
  end
  println(times[end])
end
