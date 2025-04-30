using FewBodyPhysics
using LinearAlgebra
using Plots

masses = [1e15, 1.0]
psys = ParticleSystem(masses)

K = Diagonal([0.0, 0.5])
K_transformed = psys.J * K * psys.J'

w_raw = [psys.U' * [1.0, -1.0]]
coeffs = [-1.0]

ops = [
    KineticEnergy(K_transformed);
    (CoulombPotential(c, w) for (c, w) in zip(coeffs, w_raw))...
]

a_vec = [1.0]

cfg = RunSimulationConfig(
    psys, ops, :quasirandom, 5, 1.0;
    E_exact = -0.125,
    make_gaussian = (A, w_raw) -> Rank1Gaussian(A, [a_vec])
)

run_simulation(cfg; plot_results = true)