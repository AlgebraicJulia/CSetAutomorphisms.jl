module CSetAutomorphisms

export common, canonical_hash, apply_automorphism, autos, CDict, to_vizstate,
       color_saturate, canonical_hash_nauty, init_graphs, test_iso,
       graph_to_cset

include("./Perms.jl")
include("./ColorRefine.jl")
include("./NautyInterface.jl")
include("./Canonical.jl")
include("./Visualization.jl")
include("./TestHelp.jl")

end # module