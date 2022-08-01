using Documenter
using CSetAutomorphisms

# Set Literate.jl config if not being compiled on recognized service.
config = Dict{String,String}()
if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
  config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/CSetAutomorphisms.jl/blob/gh-pages/dev"
  config["repo_root_url"] = "https://github.com/AlgebraicJulia/CSetAutomorphisms.jl/blob/main/docs"
end

makedocs(
    sitename = "CSetAutomorphisms",
    format = Documenter.HTML(),
    modules = [CSetAutomorphisms]
)


@info "Deploying docs"
deploydocs(
  target = "build",
  repo   = "github.com/AlgebraicJulia/CSetAutomorphisms.jl.git",
  branch = "gh-pages",
  devbranch = "main"
)