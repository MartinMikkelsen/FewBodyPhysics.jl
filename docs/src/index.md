# FewBodyPhysics.jl

## Installation

Get the latest stable release with Julia's package manager:

```
julia ] add FewBodyPhysics
```

## Quick example


```@example EulerSpiral
using Plots, FewBodyPhysics

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]
masses = [1, 1, 1]
K = [0 0 0; 0 1/2 0; 0 0 1/2]

```