module TestNautyInterface

using Test
using Catlab.CategoricalAlgebra, Catlab.Graphs
using CSetAutomorphisms
using CSetAutomorphisms.NautyInterface: all_autos


# Testing nauty interface
#########################
G = @acset Graph begin V=5; E=7; src=[1,2,3,2,3,1,4]; tgt=[2,3,1,1,2,3,5] end
H = @acset Graph begin V=5; E=7; src=[2,3,1,3,1,2,5]; tgt=[3,1,2,2,3,1,4] end
@test is_isomorphic(G,H)
cG, cH = call_nauty.([G,H])
@test is_isomorphic(G,cG.cset)
@test cG.cset == cH.cset
@test length(all_autos(G, cG.gens)) == cG.ngrp # this just happens to work

# Problems remain with all_autos
# sqr = @acset Graph begin V=4;E=8;src=[1,2,3,4,1,2,3,4]; tgt=[2,3,4,1,4,3,2,1] end
# cSqr = call_nauty(sqr);
# A square has D₄ symmetry id, s, r, r², r³, sr, sr², sr³
#@test cSqr.ngrp == length(all_autos(cSqr)) == 8
# Todo: seems to give incorrect result when we start with "sqr"

end # module