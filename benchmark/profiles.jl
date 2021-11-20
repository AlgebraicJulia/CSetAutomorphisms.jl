using CSetAutomorphisms, Profile, PProf
using Catlab.Graphs, Catlab.CategoricalAlgebra

# Graphs
#-------
g1 = star_graph(Graph, 5)
g2 = path_graph(Graph, 5)
g = Graph(1)
[copy_parts!(g, h) for h in [g1, g2, g1, g2]]
color_saturate(g)
@profile color_saturate(g)
pprof()
PProf.refresh(file="profile.pb.gz")