using Revise
using FewBodyPhysics
using LinearAlgebra
using Plots

masses = [1.0, 1.0, 1.0]  # proton effectively fixed
psys = ParticleSystem(masses)

K = Diagonal([1/2, 1/2, 1/2])
K_transformed = psys.J * K * psys.J'

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]

w = [psys.U' * w for w in w_list]

let
    n_basis = 200
    b1 = 6.0
    method = :quasirandom  # :quasirandom, :psudorandom
    basis_fns = GaussianBase[]
    E_trace = Float64[]
    E_best = 1e10

    for i in 1:n_basis

        bij = generate_bij(method, i, length(w_raw), b1)
        
        A = generate_A_matrix(bij, w_raw)
        push!(basis_fns, Rank0Gaussian(A))

        basis = BasisSet(basis_fns)
        ops = [KineticEnergy(K_transformed)] ∪ [CoulombPotential((w)) for w in w_raw]

        H = build_hamiltonian_matrix(basis, ops)
        S = build_overlap_matrix(basis)
        vals, _ = solve_generalized_eigenproblem(H, S)
        E0 = minimum(vals)

        push!(E_trace, E0)
        E_best = min(E_best, E0)
        println("Step $i: E = $E0")
    end

    Theortical_value = -0.2620050702328

    plot(1:n_basis, E_trace, xlabel="Number of Gaussians", ylabel="Energy [Hartree]",
         lw=2, label="SVM energy", title="Hydrogen Anion Convergence")
end
