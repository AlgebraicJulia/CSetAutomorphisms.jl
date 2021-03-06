using Catlab.CategoricalAlgebra, Catlab.Present, Catlab.Graphs
using Nauty
import LightGraphs
const lg = LightGraphs
using DataStructures: OrderedSet


"""
The general `to_lg` function is overkill when we are just trying to compute the
canonical label of a (simple) DiGraph. This output can be used with Nauty.jl by
giving an initial partition of all vertices in the same set.
"""
function graph_to_lg(g::Graph)::lg.DiGraph
  res = lg.DiGraph(nv(g))
  for (s, t) in zip(g[:src], g[:tgt])
    lg.add_edge!(res, s, t)
  end
  return lg
end

"""
Produce a simple digraph and initial coloring that encodes an ACSet. I.e. a
faithful functor from ACSet -> SimpleColoredDiGraph (for any schema C)

The fact that different vertices come from different ACSet objects can be
encoded as an initial coloring.

How do we denote different morphisms between the same pair of objects?
- subdivide the graph, putting a vertex in the middle of every edge
- color the 'edge vertex' by which morphism it corresponds to.

For example, the multi-digraph V₁ ⟶ V₂ is translated to the vertex-labeled
digraph: V₁ <- src -- E₁ -- tgt -> V₂, where the V's have the same color, all
src nodes have the same color, etc.
"""
function to_lg(g::StructACSet, p::Presentation)
  colors, curr, oinds = Any[], 1, Dict{Symbol}{UnitRange}()
  lgr = lg.DiGraph{Int}(sum(g.obs))

  # Handle Ob
  for o in [g.args[1] for g in p.generators[:Ob]]
    append!(colors, repeat([o], nparts(g, o)))
    oinds[o] = curr:length(colors)
    curr = length(colors) + 1
  end

  # Handle Hom
  for h in p.generators[:Hom]
    hom_name, d_, cd_  = h.args
    domcodom = [d_.args[1], cd_.args[1]]
    for st in enumerate(g[hom_name])
      s_ind, t_ind = [oinds[name][val] for (name, val) in zip(domcodom, st)]
      lg.add_vertex!(lgr)
      push!(colors, hom_name)
      lg.add_edge!(lgr, s_ind, lg.nv(lgr))
      lg.add_edge!(lgr, lg.nv(lgr), t_ind)
    end
    oinds[hom_name] = curr:length(colors)
    curr = length(colors) + 1
  end

  # Handle Data
  attrdict = Dict{Pair{Symbol, Any},Int}() # store the vertex for each value
  attrvals = Dict{Symbol, Vector{Any}}()   # store all values in canonical order
  for dt in p.generators[:AttrType]
    d_name, vals = dt.args[1], Set()
    # Find out how many unique elements there are of this type, across all attrs
    for h in p.generators[:Attr]
      if h.args[3].args[1] == d_name
        union!(vals, g[h.args[1]])
      end
    end
    # Give each val a unique color and a unique vertex
    attrvals[d_name] = sort(collect(vals))
    for v in attrvals[d_name]
      lg.add_vertex!(lgr)
      col = d_name => v
      push!(colors, col)
      attrdict[col] = lg.nv(lgr)
    end
    oinds[d_name] = curr:length(colors)
    curr = length(colors) + 1
  end

  # Handle Attr
  attr_colors = length(colors)+1
  for h in p.generators[:Attr]
    attr_name, d_, cd_  = h.args
    dom_,codom_ = [d_.args[1], cd_.args[1]]
    for (s, v) in enumerate(g[attr_name])
      s_ind = oinds[dom_][s]
      t_ind = attrdict[codom_ => string(v)]
      lg.add_vertex!(lgr)
      push!(colors, attr_name)
      lg.add_edge!(lgr, s_ind, lg.nv(lgr))
      lg.add_edge!(lgr, lg.nv(lgr), t_ind)
    end
    oinds[attr_name] = curr:length(colors)
    curr = length(colors) + 1
  end

  # Compute labels and partition from color data
  #---------------------------------------------
  # Get an array of arrays of nodes which are all the same color
  colorsarray = Vector{Int32}[findall(==(c), colors)
      for c in OrderedSet(colors)]
  # Assign numeric values corresponding to the colors, starting from '0'
  labels = vcat(((y -> map(x->x-1, y)).(colorsarray))...)
  # Encode that a partition ends with a '0'
  partition = vcat([begin z[end]=0; z end
                    for z in ones.(Cint,size.(colorsarray))]...)
  return (lgr, labels, partition, oinds, attrvals, colorsarray, attr_colors)
end

"""
Get canonical ACSet from Nauty's canong output
"""
function from_canong(g::Vector{UInt64}, p::Presentation, G::StructACSet
                    )::StructACSet
  _, _, _, oinds, attrdict = to_lg(G,p)
  gadj = Nauty.label_to_adj(g)

  # Each hom with the offset of its domain and codomain
  homdata = [let (hn, d, cd)=h.args;
             (hn, [oinds[x.args[1]].start-1 for x in [d,cd]]) end
             for h in p.generators[:Hom]]
  # The ordered data for each hom
  homvals = [last.(sort([[findfirst(>(0), v) - offset for (offset, v)
          in zip(offsets, [gadj[:,i], gadj[i,:]])]
          for i in oinds[h]])) for (h,offsets) in homdata]
  # Each attr with the offset of its domain and codomain
  attrdata = [let (hn, d, cd)=h.args;
             (hn, cd.args[1], oinds[x.args[1]][1]-1 for x in [d,cd]) end
             for h in p.generators[:Attr]]
  # The ordered data for each attr
  attrvals = [last.(sort([[findfirst(>(0), v) - offset for (offset, v)
              in zip(offsets, [gadj[:,i], gadj[i,:]])]
              for i in oinds[h]])) for (h,_,offsets) in attrdata]
  I = deepcopy(G)
  for (h, vs) in zip(first.(homdata), homvals)
    set_subpart!(I, h, vs)
  end
  for ((h,cd), vs) in zip(attrdata, attrvals)
    set_subpart!(I, h, [attrdict[cd][v] for v in vs])
  end
  return I
end

"""Convert to Nauty.jl input, then interpret result back into Catlab language"""
function canonical_iso_nauty(g::StructACSet, p::Presentation)::StructACSet
  lgrph, lab, prt, _, _ = to_lg(g, p)
  # println("Calling nauty with graph of order $(lg.nv(lgrph))")
  if false #lg.nv(lgrph) > 64
    # println("slow canonical iso")
    res =  canonical_iso(g)
    # println("done")
    return res
  else
    cf = coloured_digraph_canonical_form(lgrph, lab, prt)
    return from_canong(cf, p, g)
  end
end

"""Digraph isomorphism testing in Nauty.jl"""
function coloured_digraph_canonical_form(g, labels, partition)
  options = Nauty.doptionblk_mutable()
  options.defaultptn = false
  options.getcanon = true
  return Nauty.densenauty(Nauty.lg_to_nauty(g), options, labels, partition).canong
end



function dreadnaut_input(g::StructACSet)
  p = Presentation(g)
  lgr, _, _, _, _, colorsarray, ac = to_lg(g, p)
  m = lg.adjacency_matrix(lgr)
  join(["n=$(lg.nv(lgr)) g",
        join(map(1:size(m)[1]) do r
          join(string.((x->x-1).(findall(==(1),m[r,:]))) ," ")
        end, ";"),
        ". f = [$(join([join(c.-1, ",") for c in colorsarray],"|"))]",
        "c d x b z"], " ")
end

function call_nauty(g::StructACSet)
  f = tempname()
  write(f, dreadnaut_input(g))
  io,errio = IOBuffer(), IOBuffer()
  cmd = pipeline(`cat $f`, `dreadnaut`)
  cmd2 = pipeline(cmd; stdout=io, stderr=errio)
  run(cmd2)
  str = String(take!(io))
  if isempty(str)
    estr = String(take!(errio))
    #show(stdout, "text/plain",g)
    error(estr)
  end
  _ = [Base.parse(Int, x) for x in split(split(
    str[findfirst("seconds", str)[1]:end],"\n")[2], " ")[2:end]]
  # Parsing the canonical graph
  # canong = map(split(str[findfirst("0 :", str)[1]:end],"\n")) do l
  #   es = last(split(l, ":"))[1:end-1]
  #   [parse(Int, x) for x in filter(x->!isempty(x), split(es, " "))]
  # end
  split(str, "\n")[end-1]
end
