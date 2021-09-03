using Catlab.Present, Catlab.Theories, Catlab.CategoricalAlgebra

@present TheoryRxn(FreeSchema) begin
  (Molecule, Atom, Bond)::Ob
  (b₁, b₂)::Hom(Bond,Atom)
  mol_in::Hom(Atom, Molecule)

  (Float, Num, YesNo)::AttrType
  atomic_number::Attr(Atom, Num)
  coefficient::Attr(Molecule, Float)
  is_reactant::Attr(Molecule, YesNo)
end

@acset_type RxnGeneric(TheoryRxn)
Rxn = RxnGeneric{Float64, Int, Bool}
H2 = Rxn()
add_part!(H2, :Molecule, coefficient=2.0, is_reactant=false)
add_parts!(H2, :Atom, 2, atomic_number=[1,1], mol_in=[1,1])
add_part!(H2, :Bond, b₁=1, b₂=2)
O2 = deepcopy(H2)
set_subpart!(O2, :atomic_number, [8,8])
set_subpart!(O2, :coefficient, 1.0)
H2O = deepcopy(H2)
set_subpart!(H2O, :is_reactant, true)
add_part!(H2O, :Atom, mol_in=1, atomic_number=8)
add_parts!(H2O, :Bond, 2, b₁=[3,3], b₂=[1,2])
r1, r2 = Rxn(), Rxn()
[copy_parts!(r1, x) for x in [H2, O2, H2O]]
[copy_parts!(r2, x) for x in [H2O, H2, O2]]

println(r1==r2) # false
println(is_isomorphic(r1, r2)) # true


