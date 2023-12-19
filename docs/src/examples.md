# Examples 
    
Suppose you want to calculate the ground state energy of the hydrogen anion in the rest-frame of the proton. 

```@example example1
using Plots, FewBodyPhysics

masses = [Inf, 1, 1]
w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]

K = [0 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(masses)
K_transformedformed = J * K * J'
w_transformedformedformed = [U' * w_list[i] for i in 1:length(w_list)]
Theortical_value = -0.527751016523

p, Energy, bases = run_simulation(50, :quasirandom, w_transformedformedformed, K_transformedformed)
plot(p)
```

```@example example1
@show Theortical_value = -0.527751016523 - Energy
```