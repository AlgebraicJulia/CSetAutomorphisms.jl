module CSetAutomorphisms

export common, canonical_hash, apply_automorphism, autos, CDict, to_vizstate,
       color_saturate

include("./Perms.jl")
include("./ColorRefine.jl")
include("./Canonical.jl")
include("./Visualization.jl")

end # module