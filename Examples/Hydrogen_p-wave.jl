using FewBodyPhysics
using LinearAlgebra
using Plots

masses = [1e10, 1.0]
psys = ParticleSystem(masses)

K = Diagonal([0.0, 0.5])
K_transformed = psys.J * K * psys.J'

w_raw = [psys.U' * [1, -1]]
coeffs = [-1.0]

n_basis = 5
method = :quasirandom
b1 = default_b0(psys.scale)

basis_fns = GaussianBase[]
E₀_list = Float64[]

a_vec = [1.0]  

for i in 1:n_basis
    bij = generate_bij(method, i, length(w_raw), b1)
    A = generate_A_matrix(bij, w_raw)
    push!(basis_fns, Rank1Gaussian(A, [a_vec]))

    basis = BasisSet(basis_fns)
    ops = Operator[
        KineticEnergy(K_transformed);
        (CoulombPotential(c, w) for (c, w) in zip(coeffs, w_raw))...
    ]

    H = build_hamiltonian_matrix(basis, ops)
    S = build_overlap_matrix(basis)

    vals, _ = solve_generalized_eigenproblem(H, S)
    valid = vals .> 1e-8
    S⁻¹₂ = Diagonal(1 ./ sqrt.(vals[valid]))
    H̃ = S⁻¹₂ * H[valid, valid] * S⁻¹₂
    eigvals = eigen(H̃).values

    E₀ = minimum(eigvals)
    push!(E₀_list, E₀)
    println("Step $i: E₀ = $E₀")
end

E_exact = -0.125  #
E_min = minimum(E₀_list)
@show ΔE = abs(E_min - E_exact)

plot(1:n_basis, E₀_list, xlabel="Number of Gaussians", ylabel="E₀ [Hartree]",
     lw=2, label="E₀ estimate", title="p-wave Hydrogen Convergence")
hline!([E_exact], label="Exact: -0.125", linestyle=:dash)
