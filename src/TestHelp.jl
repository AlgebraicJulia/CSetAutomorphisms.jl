"""
Functions that help the automation of tests
"""

using Test
using Catlab.CategoricalAlgebra, Catlab.Present
using Catlab.CategoricalAlgebra.CSetDataStructures: struct_acset

xs(x::Int)::Symbol = Symbol("x$x")
xs(xx::AbstractVector{Int})::Vector{Symbol} = [Symbol("x$x") for x in xx]
es(x::Int)::Symbol = Symbol("e$x")
es(xx::AbstractVector{Int})::Vector{Symbol} = [Symbol("e$x") for x in xx]

"""
Create a CSet type specified by a graph
Vertices are x₁,x₂,..., edges are e₁, e₂,...
all edges are indexed
"""
function graph_to_cset(grph::StructACSet, name::Symbol)::Pair{Presentation,StructACSet}
  pres = Presentation(FreeSchema)
  xobs = [Ob(FreeSchema,xs(i)) for i in 1:nv(grph)]
  for x in xobs
    add_generator!(pres, x)
  end
  for (i,(src, tgt)) in enumerate(zip(grph[:src], grph[:tgt]))
    add_generator!(pres, Hom(es(i), xobs[src], xobs[tgt]))
  end

  expr = struct_acset(name, StructACSet, pres, index=es(1:ne(grph)))
  eval(expr)
  csettype = eval(name)
  return pres => Base.invokelatest(csettype)
end

  """Create n copies of a CSet based on a graph schema"""
function init_graphs(name::Symbol, schema::StructACSet, consts::Vector{Int},
             n::Int=2)::Pair{Presentation, Vector{StructACSet}}
  pres, cset = graph_to_cset(schema, name)
  for (i, con) in enumerate(consts)
    add_parts!(cset, Symbol("x$i"), con)
  end
  return pres => [deepcopy(cset) for _ in 1:n]
end

  """Confirm canonical hash tracks with whether two ACSets are iso"""
function test_iso(a::StructACSet,b::StructACSet, pres::Union{Nothing,Presentation}=nothing)::Test.Pass
  @test a != b  # confirm they're not literally equal
  eq = is_isomorphic(a,b)

  tst = a -> eq ? a : !a
  if !isnothing(pres)
    @test tst(canonical_hash(a,pres=pres) == canonical_hash(b; pres=pres))
  end
  @test tst(canonical_hash(a) == canonical_hash(b))
end