module CSetAutomorphisms

export common, canonical_hash, apply_automorphism, autos, CDict, to_vizstate

include("./Perms.jl")
include("./ColorRefine.jl")
include("./Canonical.jl")
include("./Visualization.jl")

end # module