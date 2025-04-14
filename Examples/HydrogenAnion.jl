using Revise
using FewBodyPhysics
using LinearAlgebra
using Plots

# --- Physical setup: H⁻ system (proton + 2 electrons) ---
masses = [1e15, 1.0, 1.0]  # proton effectively fixed
psys = ParticleSystem(masses)

# Kinetic energy matrix in Jacobi coordinates
K = Diagonal([0.0, 1/2, 1/2])
K_transformed = psys.J * K * psys.J'

# Weight vectors for Coulomb interactions
w_list = [
    [1.0, -1.0, 0.0],   # e⁻–e⁻
    [1.0,  0.0, -1.0],  # p–e⁻
    [0.0,  1.0, -1.0],  # p–e⁻
]
w_transformed = [(psys.U' * w) for w in w_list]


# --- SVM energy loop ---
let
    n_basis = 250
    b1 = 50.0
    n_terms = length(w_transformed)

    basis_fns = GaussianBase[]
    E_trace = Float64[]
    E_best = 1e10

    for i in 1:n_basis
        bij = -log.(halton(i, n_terms)) * b1
        A = generate_A_matrix(bij, w_transformed)
        new_fn = Rank0Gaussian(A)
        push!(basis_fns, new_fn)

        basis = BasisSet(basis_fns)
        ops = [KineticEnergy(K_transformed)] ∪ [CoulombPotential(w) for w in w_transformed]

        H = build_hamiltonian_matrix(basis, ops)
        S = build_overlap_matrix(basis)
        vals, _ = solve_generalized_eigenproblem(H, S)
        E0 = minimum(vals)

        push!(E_trace, E0)
        E_best = min(E_best, E0)
        println("Step $i: E = $E0")
    end

    # Report and plot
    Theoretical_value = -0.527751016523
    @show Theoretical_value - E_best

    plot(1:n_basis, E_trace, xlabel="Number of Gaussians", ylabel="Energy [Hartree]",
         lw=2, label="SVM energy", title="Hydrogen Anion Convergence")
end
