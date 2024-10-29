using Revise
using .FewBodyPhysics

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]
masses = [1.0,1.0,1.0]
K = [1/2 0 0; 0 1/2 0; 0 0 1/2]
J, U = jacobi_transform(masses)
K_transformed = J * K * J'
w_transformed = [U' * w for w in w_list]

Theortical_value = -0.2620050702328

p, Energy = run_simulation(50,:quasirandom, w_transformed, K_transformed)
@show Energy-Theortical_value
display(p)