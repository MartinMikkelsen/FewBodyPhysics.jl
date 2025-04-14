module Sampling

using ..Types
using ..Coordinates
using ..Hamiltonian

using LinearAlgebra

export generate_basis, compute_ground_state_energy

"""
generate_basis(widths::Vector{Matrix{Float64}}, rank::Int=0)

Construct a `BasisSet` from a list of correlation matrices and optional rank.
"""
function generate_basis(widths::Vector{Matrix{Float64}}, rank::Int=0)
    if rank == 0
        funcs = [Rank0Gaussian(A) for A in widths]
    else
        error("Only Rank0Gaussian implemented in generate_basis")
    end
    return BasisSet(funcs)
end

"""
compute_ground_state_energy(basis::BasisSet, ops::Vector{Operator})

Construct the Hamiltonian and overlap matrices and return the lowest eigenvalue.
"""
function compute_ground_state_energy(basis::BasisSet, ops::Vector{Operator})
    H = build_hamiltonian_matrix(basis, ops)
    S = build_overlap_matrix(basis)
    vals, _ = solve_generalized_eigenproblem(H, S)
    return minimum(vals)
end

end # module