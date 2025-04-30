module Optimization

using LinearAlgebra, Plots
using ..Types
using ..Hamiltonian
using ..MatrixElements
using ..Coordinates
using ..Sampling
export RunSimulationConfig, run_simulation

struct RunSimulationConfig
    psys::ParticleSystem
    ops::Vector{Operator}
    method::Symbol
    n_basis::Int
    b1::Float64
    E_exact::Union{Float64, Nothing}
    make_gaussian::Function
end

function RunSimulationConfig(psys, ops, method, n_basis, b1;
                             E_exact=nothing,
                             make_gaussian=(A, w_raw) -> Rank0Gaussian(A))
    return RunSimulationConfig(psys, ops, method, n_basis, b1, E_exact, make_gaussian)
end

function run_simulation(cfg::RunSimulationConfig; plot_results::Bool = true)
    basis_fns = GaussianBase[]
    E₀_list = Float64[]
    w_raw = [op.w for op in cfg.ops if op isa CoulombPotential]

    for i in 1:cfg.n_basis
        bij = generate_bij(cfg.method, i, length(w_raw), cfg.b1)
        A = generate_A_matrix(bij, w_raw)
        push!(basis_fns, cfg.make_gaussian(A, w_raw))

        basis = BasisSet(basis_fns)
        H = build_hamiltonian_matrix(basis, cfg.ops)
        S = build_overlap_matrix(basis)

        λs, Us = eigen(Symmetric(S))
        keep = λs .> 1e-12
        S⁻¹₂ = Us[:, keep] * Diagonal(1 ./ sqrt.(λs[keep])) * Us[:, keep]'
        H̃ = Symmetric(S⁻¹₂ * H * S⁻¹₂)
        vals, vecs = eigen(H̃)
        E₀ = minimum(vals)
        c₀ = vecs[:, argmin(vals)]

        push!(E₀_list, E₀)
    end

    E_min = minimum(E₀_list)
    if cfg.E_exact !== nothing
        ΔE = abs(E_min - cfg.E_exact)
        @show ΔE
    end

    if plot_results
        plot(1:cfg.n_basis, E₀_list, xlabel="Number of Gaussians", ylabel="E₀ [Hartree]",
             lw=2, label="E₀ estimate", title="Energy Convergence")
        if cfg.E_exact !== nothing
            hline!([cfg.E_exact], label="Exact", linestyle=:dash)
        end
        display(current())
    end

    return E₀_list
end



end # module
