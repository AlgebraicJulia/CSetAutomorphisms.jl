using Test

@testset "Canonical" begin
  include("Canonical.jl")
end

@testset "NautyInterface" begin
  include("NautyInterface.jl")
end

@testset "Diagrams" begin
  include("Diagrams.jl")
end
