using Documenter, FewBodyPhysics

push!(LOAD_PATH,"../src/")
makedocs(
    source  = "src", 
    #workdir = "build", 
    sitename = "FewBodyPhysics.jl",
    #format = Documenter.HTML(),)
deploydocs(
    repo = "github.com/MartinMikkelsen/FewBodyPhysics.jl.git",
    target = "docs/build",
    branch="gh-pages",
)