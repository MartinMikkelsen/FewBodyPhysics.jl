using FewBodyPhysics

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]
masses = [Inf, 1, 1]
K = [0 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(masses)
K_transformed = J * K * J'
w_transformed = [U' * w_list[i] for i in 1:length(w_list)]

Theortical_value = -0.527751016523

p, Energy = run_simulation(50, :quasirandom, w_transformed, K_transformed)

OptimizeGlobalParameters(5,1,10,masses,0)
@show Energy-Theortical_value

display(p)

