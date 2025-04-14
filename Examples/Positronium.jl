using Revise
using FewBodyPhysics
using LinearAlgebra
using Plots

# Define system
masses = [1.0, 1.0, 1.0]
psys = ParticleSystem(masses)

# Define kinetic matrix K (unitless, internal coordinates)
K = Matrix(I, 3, 3) ./ 2
K_transformed = psys.J * K * psys.J'

# Weight vectors (relative positions in particle space)
w_list = [
    [1.0, -1.0, 0.0],
    [1.0, 0.0, -1.0],
    [0.0, 1.0, -1.0],
]
w_transformed = [normalize(psys.U' * w) for w in w_list]

# Generate basis set
n_basis = 500
n_terms = length(w_transformed)  # should be 3
bij = 0.5 .+ rand(n_terms)

A = generate_A_matrix(bij, w_transformed)
basis = BasisSet([
    Rank0Gaussian(generate_A_matrix(rand(length(w_transformed)), w_transformed))
    for _ in 1:n_basis
])

# Define operators
ops = [KineticEnergy(), CoulombPotential(w_transformed[1]),
       CoulombPotential(w_transformed[2]), CoulombPotential(w_transformed[3])]

# Solve
H = build_hamiltonian_matrix(basis, ops)
S = build_overlap_matrix(basis)
vals, _ = solve_generalized_eigenproblem(H, S)
E0 = minimum(vals)

# Compare to known result
Theoretical_value = -0.2620050702328
@show E0 - Theoretical_value
