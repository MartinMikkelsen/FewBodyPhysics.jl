# FewBodyPhysics.jl

## Installation

Get the latest stable release with Julia's package manager:

```
julia ] add FewBodyPhysics
```

## Quick example

```@example
using Plots, FewBodyPhysics

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]
masses = [1,1,1]
K = [1/2 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(masses)
K_trans = J * K * J'
w_trans = [U' * w_list[i] for i in 1:length(w_list)]

Theortical_value = -0.2620050702328

p, Energy = run_simulation(50,:quasirandom, w_trans, K_trans)
plot(p) # hide
```