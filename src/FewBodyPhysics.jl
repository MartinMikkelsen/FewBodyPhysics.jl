module FewBodyPhysics

include("types.jl")
include("coordinates.jl")
include("constants.jl") 
include("matrix_elements.jl")
include("hamiltonian.jl")
include("sampling.jl")
include("optimization.jl")

using .Types
using .Coordinates
using .MatrixElements
using .Hamiltonian
using .Sampling
using .Optimization
using .Constants  

export Particle, GaussianBase, Rank0Gaussian, Rank1Gaussian, Rank2Gaussian,
       BasisSet, Operator, KineticEnergy, CoulombPotential,
       FewBodyHamiltonian, MatrixElementResult, SystemCoordinates, ParticleSystem, generate_A_matrix, run_simulation

export compute_matrix_element, build_overlap_matrix, build_operator_matrix,
       build_hamiltonian_matrix, solve_generalized_eigenproblem,
       generate_basis, compute_ground_state_energy,
       corput, halton, optimize_ground_state_energy, jacobi_kinetic_matrix

export m_p, m_n, m_pi0, m_pi, ħc, μ, m_bare 

end
