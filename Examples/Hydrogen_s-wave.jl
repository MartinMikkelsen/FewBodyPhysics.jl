using FewBodyPhysics
using LinearAlgebra
using Plots

masses = [1e15, 1.0]
psys = ParticleSystem(masses)

K = Diagonal([0.0, 1/2])
K_transformed = psys.J * K * psys.J'

w = psys.U' * [1.0, -1.0]
Operators = [
    KineticEnergy(K_transformed),
    CoulombPotential(-1.0, w)
]

FewBodySystem = RunSimulationConfig(psys, Operators, :quasirandom, 10, 1.5; E_exact = -0.5)

run_simulation(FewBodySystem)
