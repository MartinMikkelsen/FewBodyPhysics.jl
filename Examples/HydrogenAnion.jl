using Revise
using .FewBodyPhysics

masses = [1e15, 1.0, 1.0]
w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]

K = [0 0 0; 0 1/2 0; 0 0 1/2]
J, U = jacobi_transform(masses)
K_transformed = J * K * J'
w_transformed = [U' * w for w in w_list]
Theortical_value = -0.527751016523

p, Energy, bases = run_simulation(50, :quasirandom, w_transformed, K_transformed)

