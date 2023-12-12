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
masses = [1, 1, 1]
K = [0 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(masses)
K_transformed = J * K * J'
w_transformed = [U' * w_list[i] for i in 1:length(w_list)];

p = run_simulation(50, :quasirandom, w_transformed, K_transformed)
```