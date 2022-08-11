using Catlab.CategoricalAlgebra, Catlab.Present, Catlab.Graphs
using DataStructures: OrderedSet, DefaultDict


bashit(str) = run(`bash -c "$str"`)


struct CPerm
  obs::Dict{Symbol, FinFunction}
  homs::Dict{Symbol, Pair{FinFunction, FinFunction}}
end


function CPerm2(oinds::Dict{Symbol,UnitRange},canon::Vector{Int}, S)
  canon2 = Dict([k=>canon[v].-(v.start-2) for (k,v) in oinds])
  omap = Dict([k=>FinFunction(v) for (k,v) in collect(canon2) if k ∈ ob(S)])
  hmap = Dict(map(zip(hom(S),dom(S),codom(S))) do (h,s,t)
          σt,σsi = FinFunction.([canon2[t],(canon2[s])])
          h=>(σsi=>σt)
        end)
  CPerm(omap, hmap)
end

function CPerm(oinds::Dict{Symbol,UnitRange},canon::Vector{Int}, S)
  canon2 = Dict([k=>canon[v].-(v.start-2) for (k,v) in oinds])
  omap = Dict([k=>FinFunction(v) for (k,v) in collect(canon2) if k ∈ ob(S)])
  hmap = Dict(map(zip(hom(S),dom(S),codom(S))) do (h,s,t)
          σt,σsi = FinFunction.([canon2[t],(canon2[s])])
          h=>(σsi=>σt)
        end)
  CPerm(omap, hmap)
end


"""
- orb: partitions of Cset parts into orbits
  e.g. if E#2 = E#5, then these two elements are symmetric
- ngrp: number of elements in the automorphism group
- iso from nauty input into canonical isomorph
- cset: canonical isomorph, codom of cmap
"""
struct NautyRes
  hsh::String
  orb::Dict{Symbol, Vector{Int}}
  gens::Vector{CPerm} # generating perms
  ngrp::Int
  cmap::ACSetTransformation
  cset::StructACSet
end

"""
Make shell command to dreadnaut and collect output
"""
function call_nauty(g::StructACSet{S}; check=false) where S
  m,oinds,_ = to_adj(g)

  # Run nauty and collect output as a string
  inp = dreadnaut(g)
  tmp = tempname()
  dreadpth = joinpath(@__DIR__, "../deps/nauty27r3/dreadnaut")
  cmd = "echo \"$inp\" | $dreadpth > $tmp"
  bashit(cmd)
  res = open(f->read(f, String), tmp)

  # regexes
  reg_perms = r"(([ ]\d+)+\n)([ ][ ][ ]([ ]\d+)+\n)*"
  reg_canon = r"\d+ :[ ](([ ]\d+)*);"
  reg_hash = r"\[(\w+ \w+ \w+)\]"

  # parse permutations for generators + the canonical graph
  cp,gens... = reverse(map(eachmatch(reg_perms, res[1:findlast(']', res)])) do rng
    CPerm(oinds, [parse(Int, x) for x in split(strip(rng.match), r"\s+")], S)
  end)

  # parse canonical graph
  canonm = zeros(Bool, size(m))
  for (i,r) in enumerate(first.(collect.(eachmatch(reg_canon, res))))
    if !isempty(strip(r))
      canonm[i,[parse(Int,x)+1 for x in split(strip(r)," ")]] .= true
    end
  end
  canong = from_adj(g, oinds, canonm)

  # parse other things
  hshstr = match(reg_hash, res)[1]
  orb = parse_orb(g, oinds, res)
  grpsize = parse(Int, match(r"grpsize=(\d+)", res)[1])

  # sanity check
  h = apply(g,cp)
  string(codom(h)) == string(canong) || error(
    "$(codom(h))\n\n$(string(canong))")

  return NautyRes(hshstr, orb, gens, grpsize, h, canong)
end

"""Parse / postprocess orbits from the end of dreadnaut input"""
function parse_orb(g::StructACSet{S}, oinds, res::String) where S
  orbd = DefaultDict{Symbol, Vector{Vector{Int}}}(()->Vector{Int}[])
  reg = r"(\d+)(:(\d+)\s\(\d+\))?\;" # match orbits
  orb = Dict([o => zeros(Int, nparts(g, o)) for o in ob(S)])
  for m in eachmatch(reg, res[findlast(']',res):end])
    rg =[parse(Int, x)+1 for x in [m[1], isnothing(m[3]) ? m[1] : m[3]]]
    rng = collect(rg[1]:rg[2])
    symb = only([k for (k,v) in collect(oinds) if rng ⊆ v])
    push!(orbd[symb], rng)
  end
  for o in ob(S)
    for (i,v) in enumerate(orbd[o])
      orb[o][v.-(oinds[o].start-1)] .= i
    end
  end
  return orb
end

"""
Convert C-Set to symmetricadjacency matrix
"""
function to_adj(X::StructACSet{S}) where S
  p = Presentation(S)
  colors, curr, oinds = Any[], 1, Dict{Symbol}{UnitRange}()
  for o in [g.args[1] for g in p.generators[:Ob]]
    append!(colors, fill(o, nparts(X, o)))
    oinds[o] = curr:length(colors)
    curr = length(colors) + 1
  end
  n_ob = length(colors)
  mat = zeros(Bool, (n_ob,n_ob))

  for h in p.generators[:Hom]
    hom_name, d_, cd_  = h.args
    hom_name_ = Symbol("_$hom_name")
    d,cd = [d_.args[1], cd_.args[1]]
    orig_rows = size(mat)[1]
    nd = nparts(X, d)
    mat = [mat zeros(Bool, (orig_rows, 2*nd))]
    n_col = size(mat)[2]
    hom_mat = zeros(Bool, (2*nd,size(mat)[2]))
    for (i,v) in enumerate(X[hom_name])
      hom_mat[i,[oinds[cd][v],oinds[d][i]]] .= true
      hom_mat[nd+i,[oinds[d][i],orig_rows+i]] .= true
    end
    append!(colors, vcat(fill(hom_name, nparts(X, d))))
    oinds[hom_name] = curr:length(colors)
    curr = length(colors) + 1
    append!(colors, vcat(fill(hom_name_, nparts(X, d))))
    oinds[hom_name_] = curr:length(colors)
    curr = length(colors) + 1
    mat = [mat;hom_mat]
  end
  mr, mc = size(mat)
  mat = [mat zeros(Bool, mr, mr-mc)]
  colorsarray = Vector{Int}[findall(==(c), colors) for c in OrderedSet(colors)]
  (mat .|| mat', oinds, colorsarray)  # symmetrize matrix
end

"""
Symmetric adjacency matrix to C-Set.
"""
function from_adj(X::StructACSet{S}, oinds::Dict{Symbol, UnitRange},
                  m::AbstractMatrix{Bool}) where S
  Y = deepcopy(X)
  inv_dict = Dict(vcat(map(collect(oinds)) do (k,vs)
    [v=>i for (i,v) in enumerate(vs)]
  end...))
  inv_ = [inv_dict[i] for i in 1:length(inv_dict)]
  for (h, s, t) in zip(hom(S),dom(S),codom(S))
    h_ = Symbol("_$h")
    for (i,h_i) in enumerate(oinds[h_])
      src_ind, hom_ind = findall(m[h_i,:])
      src_tgt = findall(m[hom_ind,:])
      tgt_ind_ = setdiff(src_tgt, vcat([h_i,src_ind...]))
      tgt_ind = isempty(tgt_ind_) ? src_ind : only(tgt_ind_)
      set_subpart!(Y, inv_dict[src_ind], h, inv_dict[tgt_ind])
    end
  end
  Y
end

"""
Construct input for dreadnaut to compute automorphism group generators,
canonical permutation/hash, and orbits.
"""
function dreadnaut(g::StructACSet)
  m,_, colorsarray = to_adj(g)
  join(["n=$(size(m)[1]) g",
        join(map(1:size(m)[1]) do r
          join(string.((x->x-1).(findall(==(1),m[r,:]))) ," ")
        end, ";"),
        ". f = [$(join([join(c.-1, ",") for c in colorsarray],"|"))]",
        "p -m c d x b z o"], " ")
end

"""
Apply a permutation on a CSet.
"""
function apply(X::StructACSet{S}, p::CPerm; fix=false,check=false) where S
  cd = deepcopy(X)
  for (h,s,t) in zip(hom(S),dom(S),codom(S))
    σs, σt = [FinFunction((collect(x))) for x in p.homs[h]]
    σti = FinFunction(invperm(collect(σt)))
    f = FinFunction(X,h)
    m =  compose(σs, f, (fix ? [] : [σti])...)
    set_subpart!(cd, h, m)
  end
  h = ACSetTransformation(X, cd; Dict([k=>FinFunction(invperm(collect(v))) for (k,v) in p.obs])...)
  if check && !is_natural(h) error("is not natural") end
  return h
end

all_autos(nr::NautyRes) =  all_autos(nr.cset, nr.gens)


function all_autos(X::StructACSet, gens::Vector{CPerm})
  seen = Dict(string(X)=>X)
  to_check = [X=>g for g in gens]
  while !isempty(to_check)
    M, g = pop!(to_check)
    M′ = codom(apply(M, g; fix=true))
    s = string(M′)
    if !haskey(seen, s)
      seen[s] = M′
      append!(to_check, [M′=>g for g in gens])
    end
  end
  return values(seen) |> collect
end

