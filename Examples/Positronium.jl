include("../src/coordinates.jl")
include("../src/matrix_elements.jl")
include("../src/sampling.jl")
include("../src/constants.jl")

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]
masses = [5496.918, 3670.481, 206.7686]
K = [0 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(masses)
K_trans = J * K * J'
w_trans = [U' * w_list[i] for i in 1:length(w_list)]

p = run_simulation(50,:quasirandom)
display(p)