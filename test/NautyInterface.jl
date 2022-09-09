module TestNautyInterface

using Revise
using Test
using Catlab.CategoricalAlgebra, Catlab.Graphs
using CSetAutomorphisms
using CSetAutomorphisms.NautyInterface: all_autos
using CSetAutomorphisms.TestHelp: Labeled

# Testing nauty interface
#########################
G = @acset Graph begin V=5; E=7; src=[1,2,3,2,3,1,4]; tgt=[2,3,1,1,2,3,5] end
H = @acset Graph begin V=5; E=7; src=[2,3,1,3,1,2,5]; tgt=[3,1,2,2,3,1,4] end
@test is_isomorphic(G,H)
cG, cH = call_nauty.([G,H])
@test is_isomorphic(G,cG.cset)
@test cG.cset == cH.cset
# this also checks that the # of automorphisms = the count nauty computes
@test all(h->dom(h)==codom(h), all_autos(G, cG.gens))

# Problems remain with all_autos
sqr = @acset Graph begin V=4;E=8;src=[1,2,3,4,1,2,3,4];
                                 tgt=[2,3,4,1,4,3,2,1] end

# A square has D₄ symmetry id, s, r, r², r³, sr, sr², sr³
@test all(h->dom(h)==codom(h), all_autos(call_nauty(sqr)))


# ACSets
G = @acset Labeled{String} begin
  V = 4; E = 4; src = [1,2,3,4]; tgt = [2,3,4,1]; dec = ["a","b","c","d"]
end;

H = @acset Labeled{String} begin
  V = 4; E = 4; src = [1,3,2,4]; tgt = [3,2,4,1]; dec = ["a","b","c","d"]
end;

cG, cH = call_nauty.([G,H])

@test cG.hsh == cH.hsh


G1 = @acset Labeled{String} begin
  V = 1; E = 1; src = [1]; tgt = [1]; dec = ["a"]
end;

H1 = @acset Labeled{String} begin
  V = 1; E = 1; src = [1]; tgt = [1]; dec = ["b"]
end;
cG1, cH1 = call_nauty.([G1,H1])

@test cG1.hsh != cH1.hsh

end # module