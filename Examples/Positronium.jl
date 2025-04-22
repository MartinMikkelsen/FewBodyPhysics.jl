using FewBodyPhysics
using LinearAlgebra
using Plots

masses = [1.0, 1.0, 1.0]
psys = ParticleSystem(masses)

K = Diagonal([1/2, 1/2, 1/2])
K_transformed = psys.J * K * psys.J'

w_list = [[1, -1, 0], [1, 0, -1], [0, 1, -1]]
w_raw = [psys.U' * w for w in w_list]

coeffs = [+1.0, -1.0, -1.0]

let
    n_basis = 50
    b1 = default_b0(psys.scale)
    method = :psudorandom
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
    Eᵗʰ = -0.2620050702328
    ΔE = abs(E₀ - Eᵗʰ)
    @show ΔE

    #  function
    r = range(0.01, 14.0, length=400)
    ρ_r = [rval^2 * abs2(ψ₀([rval, 0.0], c₀, basis_fns)) for rval in r]
    
    p1 = plot(r, ρ_r, xlabel="r (a.u.)", ylabel="r²|ψ₀(r)|²",
              lw=2, label="r²C(r)", title="Electron-Positron Correlation Function")
    

    # Energy convergence
    p2 = plot(1:n_basis, E₀_list, xlabel="Number of Gaussians", ylabel="E₀ [Hartree]",
    lw=2, label="Ground state energy", title="Positronium Convergence")

    plot(p1, p2, layout=(2, 1))

end

