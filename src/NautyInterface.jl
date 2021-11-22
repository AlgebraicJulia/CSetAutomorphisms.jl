using Catlab.CategoricalAlgebra, Catlab.Present
using Nauty
import LightGraphs

"""produce a simple digraph and initial coloring that encodes a CSet
Encode a CSet as a vertex-colored simple digraph

The fact that different vertices come from different CSet objects can be
encoded as an initial coloring.

How do we denote different morphisms between the same pair of objects?
Thanks to Misha Lavrov on math.stackexchange, we have a path forward:
- subdivide the graph, putting a vertex in the middle of every edge
- color the 'edge vertex' by which morphism it corresponds to.

This also eliminates the need for multigraphs.
"""
function to_lg(g::StructACSet, p::Presentation
               )::Tuple{LightGraphs.DiGraph, Vector{Int32}, Vector{Int32}}
  # Construct light graph and a distinct color for each object and hom
  #-------------------------------------------------------------------
  oinds, colors = get_oinds(g, p)

  lg = LightGraphs.DiGraph(sum(g.obs))

  for (h_color_ind, h) in enumerate(p.generators[:Hom])
    h_color = length(p.generators[:Ob]) + h_color_ind
    hom_name, d_, cd_ = h.args
    domcodom = [d_.args[1], cd_.args[1]]
    for st in enumerate(g[hom_name])
      s_ind, t_ind = [oinds[name][val] for (name, val) in zip(domcodom, st)]
      LightGraphs.add_vertex!(lg)
      LightGraphs.add_edge!(lg, s_ind, LightGraphs.nv(lg))
      LightGraphs.add_edge!(lg, LightGraphs.nv(lg), t_ind)
      push!(colors, h_color)
    end
  end

  # Compute labels and partition from color data
  #---------------------------------------------
  # Get an array of arrays of nodes which are all the same colour
  colorsarray = Vector{Int32}[findall(==(c), colors) for c in Set(colors)]
  # Nauty numbers its nodes from 0
  labels = vcat(((y -> map(x->x-1, y)).(colorsarray))...)
  # Give the last node of each colour a "label" of 0, otherwise 1, as Nauty requires
  partition = vcat([begin z[end]=0; z end for z in ones.(Cint,size.(colorsarray))]...)
  return (lg, labels, partition)
end


function get_oinds(g::StructACSet, p::Presentation)
  colors, curr, oinds = Int32[], 1, Dict{Symbol}{UnitRange}()
  obs = [g.args[1] for g in p.generators[:Ob]]
  for (i, o) in enumerate(obs)
    append!(colors, repeat([i], nparts(g, o)))
    oinds[o] = curr:length(colors)
    curr = length(colors) + 1
  end
  return oinds, colors
end


"""
Try to get an automorphism back out of the canonical labeling

This doesn't currently work (it doesn't even produce an automorphism, in
general) ... maybe we need to use the permutations of the fake nodes?
"""
function from_labels(g::StructACSet, p::Presentation, clabel::Vector{Int32}
                     )::StructACSet
  oinds, _ = get_oinds(g, p)
  cdict = CDict([o=>[x - minimum(clabel[i])+1 for x in clabel[i]]
                 for (o, i) in collect(oinds)])
  return apply_automorphism(g, cdict)
end

function canonical_hash_nauty(g::StructACSet, p::Presentation)::UInt64
  lgrph, lab, prt = to_lg(g, p)
  cf = Nauty.baked_canonical_form_color(lgrph, lab, prt).canong
  return hash(cf)
end

