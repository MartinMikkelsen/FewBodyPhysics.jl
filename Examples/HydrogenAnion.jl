using FewBodyPhysics

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]
masses = [5496.918, 3670.481, 206.7686]
K = [0 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(masses)
K_trans = J * K * J'
w_trans = [U' * w_list[i] for i in 1:length(w_list)]

Theortical_value = -0.527751016523
p, Energy = run_simulation(50, :quasirandom, w_trans, K_trans)
@show Energy-Theortical_value
display(p)

