using Catlab.CategoricalAlgebra.CSets
using PermutationGroups

# Color assigned to each elem of each compoennt
const CDict = Dict{Symbol, Vector{Int}}

# The maximum color of an empty color list is 0
max0(x::Vector{Int})::Int = isempty(x) ? 0 : maximum(x)

check_auto(x::CDict)::Bool = all(map(Base.isperm, values(x)))


"""
Construct permutation σ⁻¹ such that σσ⁻¹=id
"""
function invert_perm(x::CDict)::CDict
  return Dict([k=>Base.invperm(v) for (k, v) in collect(x)])
end

"""Compose permutations"""
function compose_perms(x::CDict, y::CDict)::CDict
  function compose_comp(xs::Vector{Int},ys::Vector{Int})::Vector{Int}
    return [ys[xs[i]] for i in eachindex(xs)]
  end
  return Dict([k => compose_comp(v1,y[k]) for (k, v1) in collect(x)])
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
