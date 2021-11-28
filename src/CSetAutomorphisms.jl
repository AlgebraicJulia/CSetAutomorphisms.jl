module CSetAutomorphisms

export common, canonical_hash, apply_automorphism, autos, CDict, to_vizstate,
       color_saturate, canonical_iso_nauty, init_graphs, test_iso, to_lg,
       graph_to_cset, graph_to_lg, Labeled, TheoryDecGraph

include("./Perms.jl")
include("./ColorRefine.jl")
include("./NautyInterface.jl")
include("./Canonical.jl")
include("./Visualization.jl")
include("./TestHelp.jl")

end # module