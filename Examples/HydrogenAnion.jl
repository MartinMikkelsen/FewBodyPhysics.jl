using Revise
using FewBodyPhysics
using LinearAlgebra
using Plots

masses = [1e15, 1.0, 1.0] 
psys = ParticleSystem(masses)

K = Diagonal([0.0, 1/2, 1/2])
K_transformed = psys.J * K * psys.J'

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]

w_raw = [psys.U' * w for w in w_list]

let
    n_basis = 50
    b1 = 7.0
    method = :quasirandom  # :quasirandom, :psudorandom
    basis_fns = GaussianBase[]
    E₀_list = Float64[]
    coeffs = [-1.0, -1.0, +1.0]

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
        vals, _ = solve_generalized_eigenproblem(H, S)
        E₀_step = minimum(vals)

        push!(E₀_list, E₀_step)
        println("Step $i: E₀ = $E₀_step")

    end

    E₀ = minimum(E₀_list)
    Eᵗʰ = -0.527751016523
    ΔE = abs(E₀ - Eᵗʰ)
    @show ΔE

    plot(1:n_basis, E₀_list, xlabel="Number of Gaussians", ylabel="E₀ [Hartree]",
         lw=2, label="Ground state energy", title="Hydrogen Anion Convergence")

end

