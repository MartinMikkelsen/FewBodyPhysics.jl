module Optimization

using Optim
using ..Types
using ..Coordinates
using ..Hamiltonian

export corput, halton, optimize_ground_state_energy

"""
corput(n::Int, b::Int=3)

Generates the nth van der Corput sequence value in base `b`.
"""
function corput(n::Int, b::Int=3)
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
optimize_ground_state_energy(init_widths::Vector{Matrix{Float64}}, ops::Vector{Operator}; max_iter=100)

Uses local optimization (Nelder-Mead) on Gaussian widths to minimize ground state energy.
"""
function optimize_ground_state_energy(init_widths::Vector{Matrix{Float64}}, ops::Vector{Operator}; max_iter::Int=100)
    vecdim = size(init_widths[1], 1)
    nfuncs = length(init_widths)
    flat_init = vcat([vec(A) for A in init_widths]...)

    function objective(x)
        widths = [reshape(x[(i-1)*vecdim^2+1:i*vecdim^2], vecdim, vecdim) for i in 1:nfuncs]
        basis = generate_basis(widths)
        return compute_ground_state_energy(basis, ops)
    end

    result = optimize(objective, flat_init, NelderMead(); iterations=max_iter)
    return Optim.minimum(result)
end

end 