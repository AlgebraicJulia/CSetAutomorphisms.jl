using Catlab.Present, Catlab.Theories, Catlab.CategoricalAlgebra

@present TheoryRxn(FreeSchema) begin
  (Molecule, Atom, Bond)::Ob
  inv::Hom(Bond,Bond)
  atom::Hom(Bond, Atom)
  mol::Hom(Atom, Molecule)

  (Float, Num)::AttrType
  atomic_number::Attr(Atom, Num)
  coefficient::Attr(Molecule, Float)

  compose(inv, inv) == id(Bond)
end

@acset_type RxnGeneric(TheoryRxn)
Rxn = RxnGeneric{Float64, Int}

H2 = Rxn()
add_part!(H2, :Molecule, coefficient=2.0)
add_parts!(H2, :Atom, 2, atomic_number=[1,1], mol=[1,1])
add_parts!(H2, :Bond, 2, atom=[1,2], inv=[2,1])

O2 = deepcopy(H2)
set_subpart!(O2, :coefficient, 1.0)
set_subpart!(O2, :atomic_number, [8,8])

H2O = deepcopy(H2)
set_subpart!(H2O, :coefficient, -2.0)
add_part!(H2O, :Atom, mol=1, atomic_number=8)
add_parts!(H2O, :Bond, 4, atom=[1,3,2,3], inv=[4,3,6,5])

r1, r2 = Rxn(), Rxn()
[copy_parts!(r1, x) for x in [H2, O2, H2O]]
[copy_parts!(r2, x) for x in [H2O, H2, O2]]

println(r1==r2) # false
println(is_isomorphic(r1, r2)) # true


