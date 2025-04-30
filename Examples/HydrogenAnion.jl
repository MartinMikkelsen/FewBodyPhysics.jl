using FewBodyPhysics
using LinearAlgebra
using Plots

masses = [1e15, 1.0, 1.0]
psys = ParticleSystem(masses)

K = Diagonal([0.0, 0.5, 0.5])
K_transformed = psys.J * K * psys.J'

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]
w_raw = [psys.U' * w for w in w_list]
coeffs = [-1.0, -1.0, +1.0]

ops = [
    KineticEnergy(K_transformed);
    (CoulombPotential(c, w) for (c, w) in zip(coeffs, w_raw))...
]

cfg = RunSimulationConfig(
    psys, ops, :quasirandom, 50, default_b0(psys.scale);
    E_exact = -0.527751016523
)

run_simulation(cfg)
