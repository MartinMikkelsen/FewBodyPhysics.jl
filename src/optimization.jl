module Optimization

using Optim
using ..Types
using ..Coordinates
using ..Hamiltonian

export optimize_ground_state_energy

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