module TestDiagrams

using Test
using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Theories
using CSetAutomorphisms


# test span homomorphisms and hashing
apX = path_graph(Graph, 2)
L = @acset Graph begin V=2;E=2;src=1;tgt=[1,2] end
R = star_graph(Graph, 4)
apY = star_graph(Graph, 3)

S1 = Span(homomorphism(apX,L), homomorphism(apX,R))
S2 = Span(homomorphism(apY,L), homomorphism(apY,R))

h = homomorphism(S1,S2)

P = @present Sch(FreeCategory) begin
   (X1,X2,X3)::Ob; f::Hom(X1,X2); g::Hom(X1,X3)
end
Grph = ACSetCat{Graph}()
D = FinCat(P)
F = FinDomFunctor(Dict(:X1=>L,:X2=>L,:X3=>L),Dict(:f=>id(L), :g=>id(L)), D, Grph)
@test is_functorial(F)

call_nauty(F)
ob_map(F,:X1)

end # module
