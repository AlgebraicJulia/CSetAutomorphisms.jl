module Perms

using Catlab.CategoricalAlgebra.CSets
using Catlab.Theories
using PermutationGroups

# Color assigned to each elem of each compoennt
const CDict = Dict{Symbol, Vector{Int}}

# The maximum color of an empty color list is 0
max0(x::Vector{Int})::Int = isempty(x) ? 0 : maximum(x)

check_auto(x::CDict)::Bool = all(Base.isperm, values(x))


"""
Construct permutation σ⁻¹ such that σσ⁻¹=id
"""
function invert_perms(x::CDict)::CDict
  return Dict([k=>Base.invperm(v) for (k, v) in collect(x)])
end

is_perms(x::CDict)::Bool = all([Base.isperm(v) for v in values(x)])

function compose_comp(xs::Vector{Int},ys::Vector{Int})::Vector{Int}
  return [ys[xs[i]] for i in eachindex(xs)]
end

"""Compose permutations"""
function compose_perms(x::CDict, y::CDict)::CDict
  Dict([k => compose_comp(v1,y[k]) for (k, v1) in collect(x)])
end


"""Convert sorted permutation into a single permutation"""
function to_perm(d::CDict, csum::Vector{Int})::Perm
  res = vcat([[vi + offset for vi in v]
              for ((_, v), offset) in zip(sort(collect(d)), csum)]...)
  return Perm(res)
end

"""Convert single permutation back into sorted permutation"""
function from_perm(p::Perm, syms::Vector{Symbol}, lens::Vector{Int},
                   csum::Vector{Int})::CDict
  res = CDict()
  for (k, l, off) in zip(syms, lens, csum)
    res[k] = [p[i+off]-off for i in 1:l]
  end
  return res
end

"""Enumerate the elements of a permutation group from its generators"""
function all_perms(perm_gens::Vector{CDict})::Set{CDict}
  syms = sort(collect(keys(perm_gens[1])))
  lens = map(length, [perm_gens[1][x] for x in syms])
  csum = vcat([0], cumsum(lens)[1:end-1])
  _,_,Cs = schreier_sims([to_perm(g, csum) for g in perm_gens])
  result = Set{CDict}()
  """
  Recursive function based on (sub)chain C and partial product r
  Algorithm from Fig II.1 of Alexander Hulpke's notes on Computational
  Group Theory (Spring 2010).
  https://www.math.colostate.edu/~hulpke/CGT/cgtnotes.pdf
  """
  function enum(i::Int, r)::Nothing
    leaf = i == length(Cs)
    C = Cs[i]
    for d ∈ C.orb
      xr = C[d]*r
      leaf ? push!(result, from_perm(xr, syms, lens, csum)) : enum(i+1, xr)
    end
  end

  enum(1, Perm(1:sum(lens)))

  return result
end

"""Apply a coloring to a C-set to get an isomorphic cset"""
function apply_automorphism(c::StructACSet{S}, d::CDict)::StructACSet{S} where {S}
  check_auto(d) || error("received coloring that is not an automorphism: $d")
  new = deepcopy(c)
  tabs, arrs, srcs, tgts = ob(S), hom(S), dom(S), codom(S)
  for (arr, src, tgt) in zip(arrs,srcs,tgts)
    set_subpart!(new, d[src], arr, d[tgt][c[arr]])
  end
  return new
end

"""
Suppose we have a set of candidates for "best permutation"
  (results in lexicographic min ACSet so far)

We get a new candidate that is a refinement of a previous candidate

Handle ambiguity by computing best case and worst case scenarios.

Returns
 - nothing: the two are not comparable
 - true: a < b
 - false: b < a
"""
const VUNI = Vector{Union{Nothing, Int}}

function compare_perms(a::Vector{VUNI}, b::Vector{VUNI})::Union{Nothing, Bool}
  for (va, vb) ∈ zip(a,b)
    for (ia, ib) in zip(va, vb)
      if isnothing(ia) || isnothing(ib)
        return nothing
      elseif ia < ib
        return true
      elseif ib < ia
        return false
      end
    end
  end
end

"""
Compute as much of the canonical data that will be ordered given a partial
coloring of an ACSet.

Given the data of a morphism f: A→B, the i'th element of f is: σᵦ(f(σₐ(i))).
With an arbitrary partition, σₐ(i) a subset of elements in |A|. Likewise for |B|
If the partitioning is not discrete, there will be cases where σₐ is not defined
We put nothing in this case.

If the best-case-scenario is worse than SOME automorphism's worst-case-scenario,
then this branch can be pruned.
"""
function order_perms(c::StructACSet, d::CDict, symorder::Vector{Symbol}
                    )::Pair{Vector{VUNI}, Vector{VUNI}}
  best, worst = map(collect, zip([order_perm(c,d,s) for s in symorder]...))
  best == worst || error("best $best worst $worst")
  return best => worst
end

function order_perm(c::StructACSet{S}, coloring::CDict, sym::Symbol
                    )::Pair{VUNI, VUNI} where {S}
  ind = findfirst(==(sym), hom(S))
  d, cd = dom(S)[ind], codom(S)[ind]
  res = [let x=perm_best_worst(c,coloring, sym, i, d, cd);
         isnothing(x) ? nothing=>nothing : x end for i in parts(c, cd)]
  bst, wrst = isempty(res) ? VUNI()=>VUNI() : map(collect, zip(res...))
  return bst => wrst
end

function perm_best_worst(c::StructACSet, coloring::CDict, f::Symbol, i::Int,
                         dm::Symbol, cdm::Symbol)::Union{Nothing, Pair{Int,Int}}
  σᵦinv = findall(==(i), coloring[cdm])
  finv = Set(incident(c, σᵦinv, f))
  σₐinv = findall(∈(finv), coloring[dm])
  return isempty(σₐinv) ? nothing : minimum(σₐinv) => maximum(σₐinv)
end

"""
Get vector listing nontrivial colors (which component and which color index) as
well as how many elements have that color. E.g. for (V=[1,1,2], E=[1,2,2,2,3,3])
we would get `[2=>(:V,1), 3=>(:E,2), 2=>(:E, 3)]`
"""
function get_colors_by_size(coloring::CDict)::Vector{Pair{Int,Tuple{Symbol, Int}}}
  res = []
  for (k, v) in coloring
    for color in 1:max0(v)
      n_c = count(==(color), v)
      if n_c > 1
        # Store which table and which color
        push!(res, n_c => (k, color))
      end
    end
  end
  return res
end

end # module