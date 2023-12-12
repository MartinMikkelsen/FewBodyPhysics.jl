using Documenter, FewBodyPhysics

push!(LOAD_PATH,"../src/")
makedocs(
    source  = "src", 
    #workdir = "build", 
    sitename = "FewBodyPhysics.jl",
    #format = Documenter.HTML(),
    #root = joinpath(dirname(pathof(FewBodyPhysics)), "..", "docs"),
)
deploydocs(
    repo = "github.com/MartinMikkelsen/FewBodyPhysics.jl.git",
    target = "build",
    branch="gh-pages",
)