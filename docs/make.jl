using Documenter, Pkg
using FewBodyPhysics


push!(LOAD_PATH,"../src/")
makedocs(
    source  = "src", 
    sitename = "FewBodyPhysics.jl",
    modules = [FewBodyPhysics], 

)
deploydocs(
    repo = "github.com/MartinMikkelsen/FewBodyPhysics.jl.git",
    target = "build",
    branch="gh-pages",
)