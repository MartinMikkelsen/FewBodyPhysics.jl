using Documenter
using FewBodyPhysics

makedocs(
    source  = "src", 
    workdir = "build", 
    sitename = "FewBodyPhysics.jl",
    format = Documenter.HTML(),
)
deploydocs(
    repo = "github.com/MartinMikkelsen/FewBodyPhysics.jl.git",
    target = "docs/build",
)