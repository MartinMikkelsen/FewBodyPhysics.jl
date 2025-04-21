using Revise
using FewBodyPhysics
using LinearAlgebra
using Plots

masses = [5496.918, 3670.481, 206.7686]  
psys = ParticleSystem(masses; scale=:nuclear)

K_cartesian = Diagonal(1 ./ masses)
K_transformed = psys.J * K_cartesian * psys.J'


w_list = [
    [1, -1, 0],  # t – d
    [1, 0, -1],  # t – μ⁻
    [0, 1, -1]   # d – μ⁻
]
w_raw = [psys.U' * w for w in w_list]
coeffs = [+1.0, +1.0, -1.0]  # [t–d, t–μ, d–μ]

let
    n_basis = 500
    b1 = 0.25
    method = :quasirandom
    basis_fns = GaussianBase[]
    E₀_list = Float64[]

    for i in 1:n_basis
        bij = generate_bij(method, i, length(w_raw), b1)
        A = generate_A_matrix(bij, w_raw)
        push!(basis_fns, Rank0Gaussian(A))

        basis = BasisSet(basis_fns)
        ops = Operator[
            KineticEnergy(K_transformed);
            (CoulombPotential(c, w) for (c, w) in zip(coeffs, w_raw))...
        ]

        H = build_hamiltonian_matrix(basis, ops)
        S = build_overlap_matrix(basis)
        vals, vecs = solve_generalized_eigenproblem(H, S)

        global c₀ = vecs[:, 1]  
        E₀ = minimum(vals)
        push!(E₀_list, E₀)

        println("Step $i: E₀ = $E₀")
    end

    E₀ = minimum(E₀_list)
    Eᵗʰ = -111.364511474  
    ΔE = abs(E₀ - Eᵗʰ)
    @show ΔE

    plot(1:n_basis, E₀_list,
        xlabel="Number of Gaussians", ylabel="E₀ [Hartree]",
        lw=2, label="Ground state energy", title="tdμ Ground State Convergence")
end
