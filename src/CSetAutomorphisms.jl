module CSetAutomorphisms
using Reexport

include(joinpath(@__DIR__, "Perms.jl"))
include(joinpath(@__DIR__, "ColorRefine.jl"))
include(joinpath(@__DIR__, "Canonical.jl"))
include(joinpath(@__DIR__, "NautyInterface.jl"))
include(joinpath(@__DIR__, "TestHelp.jl"))
include(joinpath(@__DIR__, "Visualization.jl"))
include(joinpath(@__DIR__, "Diagrams.jl"))


@reexport using .Perms
@reexport using .ColorRefine
@reexport using .Canonical
@reexport using .NautyInterface
@reexport using .TestHelp
@reexport using .Visualization
@reexport using .Diagrams

end # module
