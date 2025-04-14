using Revise
using FewBodyPhysics
using LinearAlgebra
using Plots

# Define system: fixed proton + two electrons
masses = [1e15, 1.0, 1.0]
psys = ParticleSystem(masses)

# Define kinetic matrix with fixed nucleus
K = Diagonal([0.0, 1/2, 1/2])
K_transformed = psys.J * K * psys.J'

# Weight vectors (proton–electron and e–e interactions)
w_list = [
    [1.0, -1.0, 0.0],
    [1.0, 0.0, -1.0],
    [0.0, 1.0, -1.0],
]
w_transformed = [normalize(psys.U' * w) for w in w_list]

# Generate basis set
n_basis = 250
n_terms = length(w_transformed)
bij = 0.5 .+ rand(n_terms)  

basis = BasisSet([
    Rank0Gaussian(generate_A_matrix(rand(n_terms), w_transformed))
    for _ in 1:n_basis
])

# Define Hamiltonian operators
ops = [KineticEnergy()] ∪ [CoulombPotential(w) for w in w_transformed]

# Solve eigenproblem
H = build_hamiltonian_matrix(basis, ops)
S = build_overlap_matrix(basis)
vals, _ = solve_generalized_eigenproblem(H, S)
E0 = minimum(vals)

# Compare to Hartree-limit for H⁻
Theoretical_value = -0.527751016523
@show E0 - Theoretical_value
