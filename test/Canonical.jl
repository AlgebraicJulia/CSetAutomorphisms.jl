module TestCanonical

# using Revise
using Test
using Catlab.CategoricalAlgebra, Catlab.Present, Catlab.Theories, Catlab.Graphs
using CSetAutomorphisms
using CSetAutomorphisms.NautyInterface: all_autos
using CSetAutomorphisms.Canonical: common
using CSetAutomorphisms.Perms: apply_automorphism
using CSetAutomorphisms.TestHelp: test_iso, init_graphs, Labeled


using Random


# Auxillary function tests
##########################

@test common([],[]) == 0
@test common([1],Int[]) == 0
@test common([1],Int[1]) == 1
@test common([1],Int[1,2]) == 1
@test common([1,2,3],Int[1,2]) == 2

# Helper functions for writing automorphism tests
#################################################

# Tests
#######
@present ThTri(FreeSchema) begin
  (V,T)::Ob
  (x1,x2,x3)::Hom(T,V)
end
@acset_type Tri(ThTri)

T1 = @acset Tri begin
  V=2; T=4;x1=[1,2,2,1];x2=[1,2,1,2];x3=[1,1,2,2]
end
T2 = deepcopy(T1)
set_subpart!(T2, :x1, [2,1,1,2])

test_iso(T1,T2)

canonical_hash(Graph())
canonical_hash(Graph(1))


G,H = Graph(4), Graph(4);
add_edges!(G,[1,2,4,4,3],[2,4,3,3,2]);
add_edges!(H,[2,3,1,4,4],[1,1,4,3,3]);
test_iso(G,H) # 196 automorphisms

Triangle = Graph(3) # f;g = h
add_edges!(Triangle, [1,1,2], [2,3,3]) # f,h,g

Tri_, (G,H) = init_graphs(:Triang, Triangle,[2,2,2])
for i in 1:3 set_subpart!(G, Symbol("e$i"), [1,1]) end
for i in 1:3 set_subpart!(H, Symbol("e$i"), [2,2]) end
test_iso(G, H)

Loop = Graph(1)
add_edge!(Loop, 1, 1)
Loo_, (G, H) = init_graphs(:Loo, Loop, [3])
set_subpart!(G, Symbol("e1"), [3,2,1])
set_subpart!(H, Symbol("e1"), [1,3,2])
test_iso(G, H)

cyclel, cycler = Graph(3), Graph(3)
add_edges!(cyclel,[1,2,3],[2,3,1])
add_edges!(cycler,[3,2,1],[2,1,3])
test_iso(cyclel, cycler)

n = 10
c1 = cycle_graph(Graph, 3*n)
pn = path_graph(Graph, n)
c2 = pn ⊕ pn ⊕ pn
add_edge!(c2, n, 2*n+1)
add_edge!(c2, 2*n, 1)
add_edge!(c2, 3*n, n+1)

test_iso(c1, c2)


Loop2 = Graph(1)
add_edges!(Loop2, [1,1],[1,1])

Loo2_, (G,H) = init_graphs(:Loo2, Loop2, [2])
set_subpart!(G, :e1, [2,1])
set_subpart!(G, :e2, [2,1])
set_subpart!(H, :e1, [1,1])
set_subpart!(H, :e2, [2,2])
test_iso(G, H)

# Example from Hartke and Radcliffe exposition of Nauty.
# G is their optimal ordering. H is the original.
G, H = Graph(9), Graph(9)
add_edges!(G,[1,1,2,2,3,3,4,4,5,6,7,8],
             [7,8,5,6,6,8,5,7,9,9,9,9])
add_edges!(H,[1,1,3,3,7,7,9,9,2,4,6,8],
             [2,4,2,6,4,8,6,8,5,5,5,5])
res, tree = autos(H);
# When branching is restricted to :V as is the case in Nauty
# length should be 13 without auto pruning
# length is 10 with auto pruning tactic #1
# length is 6 with auto pruning tactic #2 too
# However, we can branch on :E too. This leads to just length 4 soln.

test_iso(G,H) # nauty doesn't seem to terminate

"""Graph corresponding to schema for finite limit sketch for categories"""
catschema = @acset Graph begin
  V = 7
  E = 17
  src = [2,2,1, 3,3,3, 4,4,5,5,4, 6,6,6, 7,7,7]
  tgt = [1,1,2, 2,2,2, 2,2,2,2,5, 1,2,3, 1,2,3]
end
random_perm = Dict([:V=>randperm(7), :E=>randperm(17)])
catschema2 = apply_automorphism(catschema, random_perm)
test_iso(catschema,catschema2)


# ACSet Tests - to be supported by call_nauty in the future
###########################################################
G = @acset Labeled{String} begin
  V = 4; E = 4; src = [1,2,3,4]; tgt = [2,3,4,1]; dec = ["a","b","c","d"]
end;


H = @acset Labeled{String} begin
  V = 4; E = 4; src = [1,3,2,4]; tgt = [3,2,4,1]; dec = ["a","b","c","d"]
end;

test_iso(G,H)

I = @acset Labeled{String} begin
  V = 4; E = 4; src = [1,2,3,4]; tgt = [2,3,4,1]; dec = ["b","c","d","a"]
end;

N = @acset Labeled{String} begin
  V = 4; E = 4; src = [1,2,3,4]; tgt = [2,3,4,1]; dec = ["a","a","b","c"]
end;
test_iso(I,N)

K = @acset Labeled{String} begin
  V = 4; E = 4; src = [1,3,2,4]; tgt = [2,3,4,1]; dec = ["a","d","b","c"]
end;
test_iso(I,K)

G1 = @acset Labeled{String} begin
  V = 1; E = 1; src = [1]; tgt = [1]; dec = ["a"]
end;

H1 = @acset Labeled{String} begin
  V = 1; E = 1; src = [1]; tgt = [1]; dec = ["b"]
end;
test_iso(G1,H1)

end # module