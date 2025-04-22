module Sampling

using ..Types
using ..Coordinates
using ..Hamiltonian
using LinearAlgebra

export corput, halton, generate_basis, compute_ground_state_energy, generate_bij

"""
    corput(n::Int, b::Int=2)

Generates the nth van der Corput sequence value in base `b`. The default base is 2, which provides good uniformity for the first few dimensions and is standard in quasi-Monte Carlo settings.
"""
function corput(n::Int, b::Int=2)
    q, bk = 0.0, 1 / b
    while n > 0
        n, rem = divrem(n, b)
        q += rem * bk
        bk /= b
    end
    return q
end

"""
halton(n::Int, d::Int)

Generates the nth Halton sequence vector in `d` dimensions.
"""
function halton(n::Int, d::Int)
    base = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281]
    @assert d <= length(base) "Dimension exceeds available bases."
    return [corput(n, base[i]) for i in 1:d]
end

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
    generate_bij(method::Symbol, i::Int, n_terms::Int, b1::Float64) -> Vector{Float64}

Generate a bij vector using the specified sampling method.
Supported methods: :quasirandom, :psudorandom
"""
function generate_bij(method::Symbol, i::Int, n_terms::Int, b1::Float64)
    if method == :quasirandom
        return halton(i, n_terms) * b1
    elseif method == :psudorandom
        return rand(n_terms) * b1
    else
        error("Unsupported sampling method: $method")
    end
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