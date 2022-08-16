module Diagrams
export homomorphism, DiagCSet, ACSetCat

using Catlab.CategoricalAlgebra
import Catlab.CategoricalAlgebra: homomorphism

import ..NautyInterface: call_nauty, NautyRes

using DataStructures: DefaultDict

"""
A span homomorphism is a homomorphism between apexes and isomorphisms
between the feet of the span. This function is a work in progress.

       X
    ↙  |  ↘
  X₁   |    X₂
≅ ↓    |    ↓ ≅
  Y₁   |    Y₂
    ↖  v  ↗
       Y

TODO: this isn't implemented correctly, currently, for two reasons.

1. We need to search over all possible combinations of isomorphisms among the
   legs, as opposed to just picking an arbitrary one.

2. Each part of X maps to a part in Yₙ for each n. The preimage along Y's leg #n
   gives a set of acceptable values that this part of X can map to. The code
   currently stops if there is any preimage that has more than one element
   (because this is a constraint that has already been implemented in
   homomorphism search) but to be correct we must implement search where each
   x∈X can be designated a list of allowable values in Y to map to.

The relevance of this code being in this repo is that it seems likely that
CSetAutomorphisms can help us avoid the combinatorial explosion of searching
over all combinations of isomorphisms between the pairs of legs. It is not as
simple as converting each of the leg C-Sets into canonical form: imagine
x₁ ↦ y₁₁ and y₁↦ y₁₂ ... ignoring isomorphisms, it would seem like we cannot
map x₁ ↦ y₁. However, if y₁₁ and y₁₂ are in the same *orbit*, then there does
exist an isomorphism that would allow us to do this mapping. The challenge is
when we then continue for x₂, etc. We've actually broken some symmetry by
distinguishing y₁₁. But it seems like color saturation would allow us to break
whatever orbits are broken by this distinction and then continue in the same
way with the rest of the elements in X. This is future work.
"""
function homomorphism(X::Multispan,
                      Y::Multispan; kw...)
  getS(::StructACSet{S}) where S = S
  S = getS(apex(X)) # for some reason, type dispatch doesn't find this method
                    # if we use Multispan{StructACSet{S}}
  isos = [isomorphism(codom.(lxy)...) for  lxy in zip(legs(X), legs(Y))]
  if any(isnothing, isos) return nothing end
  init = DefaultDict(()->Dict{Int,Int}())
  for o in ob(S)
    for p in parts(apex(X),o)
      for (lx,ly,i) in zip(legs(X), legs(Y), isos)
        tgt = compose(lx[o],i[o])(p)
        preim = preimage(ly[o], tgt)
        if length(preim) > 1
          # in reality, we need to maintain a set of possible things that p can
          # map to and intersect that list with the preimage for every leg
          return nothing
        elseif length(preim) == 1
          tgt = only(preim)
          prev_tgt = get(init[o],p)
          if !isnothing(prev_tgt) && prev_tgt != tgt return nothing end
          init[o][p] = tgt
        end
      end
    end
  end
  f = homomorphism(apex(X), apex(Y); initial=init, kw...)
  if isnothing(f) return nothing end
  return f => isos # maybe define a SpanTransformation struct?
end


const ACSetCat{S} = TypeCat{S, ACSetTransformation}
const DiagCSet{D,S} = FinFunctor{D, ACSetCat{S}}

"""Invert a C-Set isomorphism"""
inv(f::CSetTransformation) = CSetTransformation(dom(f),codom(f);
  Dict(map(collect(pairs(components(f)))) do (k,v)
    k=>invperm(collect(v))
  end)...) # I think Evan has a helper function to do map pairs in one line

"""
A diagram in a C-Set category can be given a canonical form induced by the
canonical forms of the underlying C-Sets. For example, we can get a canonical
form for spans by considering diagrams of shape • ⟵ • ⟶ •
"""
function call_nauty(X::DiagCSet{D_,S}) where {D_,S}
  D = dom(X)
  # Automorphism data of all underlying objects
  canon_obs = Dict([o=>call_nauty(ob_map(X,o)) for o in ob_generators(D)])
  # Canonical form of all underlying objects
  canon_obs_ = [k=>c.cset for (k,c) in collect(canon_obs)]
  # Updating homomorphisms to point between the canonical C-Sets
  canon_hom = Dict(map(hom_generators(D)) do h
    s, t = dom(D,h), codom(D,h)
    h => compose(inv(canon_obs[s].cmap), hom_map(X,h), canon_obs[t].cmap)
  end)
  # New canonical diagram
  canon_diag = FinDomFunctor(canon_obs_, canon_hom, dom(X), codom(X))
  # orbit data is per-object
  orb = [k=>c.orb for (k,c) in collect(canon_obs)]
  # automorphism generator data is per-object
  gen = [k=>c.gens for (k,c) in collect(canon_obs)]
  # Total # of automorphisms is the product of each objects auto group
  ngrp = prod([c.ngrp for c in values(canon_obs)])
  # Construct a diagram morphism (w/ id shape map) from the input data
  dh = FinTransformation(Dict([o=>c.cmap for (o,c) in collect(canon_obs)]),
                          X, canon_diag)
  # Hash the object data and the morphism data
  new_hsh = [[l.hsh for l in values(canon_obs)]..., map(values(canon_hom)) do l
              Dict([k=>collect(v) for (k,v)
                   in collect(pairs(components(l)))]) end...]

  return NautyRes(new_hsh,orb,gen,ngrp, dh, canon_diag)
end

end # module