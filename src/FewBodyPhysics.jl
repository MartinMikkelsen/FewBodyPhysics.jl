module FewBodyPhysics

using LinearAlgebra, Plots

include("coordinates.jl")
include("matrix_elements.jl")
include("sampling.jl")

p = run_simulation(15,:psudorandom)
display(p)
end