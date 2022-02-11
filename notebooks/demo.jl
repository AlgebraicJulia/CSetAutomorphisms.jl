### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6e54a3da-863b-11ec-1cdd-c3289af2a58f
begin 
	using Revise, Pkg, PlutoUI
	Pkg.develop(path="..")
	Pkg.develop(path="../../Catlab.jl")
	using CSetAutomorphisms
	using Catlab.Graphics, Catlab.Graphs, Catlab.WiringDiagrams, Catlab.CategoricalAlgebra
	using Catlab.Present, Catlab.Theories
	using Catlab.Graphics: to_graphviz, Graphviz, TopToBottom
	md"""(imports)"""

end

# ╔═╡ 16e47677-654a-456c-8ea3-83b146f6f94b
html"""<style>
main {
    max-width: 65%;
	align: "right";
	
}
pluto-helpbox {
    display: none !important;
    width: 0px !important;
    height: 0px !important;
    padding: 0 !important;

}
"""

# ╔═╡ b9fd2d59-d95e-4e33-b87a-6797ff7da6fd
show_diagram(d::WiringDiagram, msg::String, 	 
             node_colors::Dict{Int,String})=to_graphviz(
    d,orientation=TopToBottom,labels=false,label_attr=:xlabel, 
    edge_attrs=Graphviz.Attributes(:fontname=>"Courier"),  
    title=Graphviz.Label("t", msg), node_colors=node_colors,
    node_attrs=Graphviz.Attributes(:fontname=>"Courier"))

# ╔═╡ e95f144b-4ac8-4658-bacc-b6c83c80e12c
begin 
	@present ThTri(FreeSchema) begin
  (V,T)::Ob
  (x1,x2,x3)::Hom(T,V)
end
#@acset_type Tri(ThTri)

T1 = Tri()
add_parts!(T1, :V, 2); add_parts!(T1, :T, 4)
set_subpart!(T1, :x1, [1,2,2,1])
set_subpart!(T1, :x2, [1,2,1,2])
set_subpart!(T1, :x3, [1,1,2,2])
T2 = deepcopy(T1)
set_subpart!(T2, :x1, [2,1,1,2])
end;

# ╔═╡ 6ce09a9f-b4d8-439b-bdf9-6e120ee7deb8
T2

# ╔═╡ 566ae365-767d-4a18-b7b3-8e9895a40534
isos = isomorphisms(T1,T2)

# ╔═╡ 50a94b27-d0c9-442d-8805-08ad67966e96
autos(T1)[1]

# ╔═╡ 35c4cd4d-80b7-4cb0-bd4d-a1b12cd3eea1
ag2 = collect(autos(T2)[1])

# ╔═╡ 94086b94-811b-4417-a1f9-e3fbd9005789
T2s = apply_automorphism(T2, only(ag2))

# ╔═╡ b4665060-e43c-4931-aa02-bcb26faa8615
canonical_iso(T2)

# ╔═╡ 1ae76726-31a8-4254-ad55-e2e14b8efa52
canonical_iso(T1)

# ╔═╡ 47fc37d6-9f88-40c0-b869-afbbf9873db4
canonical_hash(T2) == canonical_hash(T1)

# ╔═╡ 59003a72-9099-4da1-9409-6da1d8b9db21
hists = autos(T2;history=true)[3];

# ╔═╡ 8249440a-4a06-4583-a677-285164f26f3e
vs = to_vizstate(hists);

# ╔═╡ 377a6ace-6008-4383-9c70-904b5180b2a0
@bind n Slider(1:length(vs); show_value=true)

# ╔═╡ a2538399-f7db-4850-afd5-813c172ef03e
let v = vs[n]; show_diagram(v.wd, v.msg, v.colors) end

# ╔═╡ Cell order:
# ╠═6e54a3da-863b-11ec-1cdd-c3289af2a58f
# ╟─16e47677-654a-456c-8ea3-83b146f6f94b
# ╠═b9fd2d59-d95e-4e33-b87a-6797ff7da6fd
# ╠═e95f144b-4ac8-4658-bacc-b6c83c80e12c
# ╠═6ce09a9f-b4d8-439b-bdf9-6e120ee7deb8
# ╠═566ae365-767d-4a18-b7b3-8e9895a40534
# ╠═50a94b27-d0c9-442d-8805-08ad67966e96
# ╠═35c4cd4d-80b7-4cb0-bd4d-a1b12cd3eea1
# ╠═94086b94-811b-4417-a1f9-e3fbd9005789
# ╠═b4665060-e43c-4931-aa02-bcb26faa8615
# ╠═1ae76726-31a8-4254-ad55-e2e14b8efa52
# ╠═47fc37d6-9f88-40c0-b869-afbbf9873db4
# ╠═59003a72-9099-4da1-9409-6da1d8b9db21
# ╠═8249440a-4a06-4583-a677-285164f26f3e
# ╟─377a6ace-6008-4383-9c70-904b5180b2a0
# ╟─a2538399-f7db-4850-afd5-813c172ef03e
