# FewBodyPhysics.jl

## Installation

Get the latest stable release with Julia's package manager:

```
julia ] add FewBodyPhysics
```

## Example

We consider a system of a positron and two electrons. The energy of this system has been very accurately calculated by various approaches and it has been found to be -0.262005 in atomic units (a.u.). We calculate the ground-state energy of this systems using correlated Gaussian bases constructed stochastically with pseudorandom and quasirandom sequences. The Hamiltonian of the system is given by
```math
H = - \sum_{i=1}^{3} \frac{1}{2m_i}\frac{\partial^2}{\partial \boldsymbol{r}_i^2} + \sum_{i<j=1}^3 \frac{q_i q_j}{|\boldsymbol{r}_i-\boldsymbol{r}_j|}
```
The masses of the three constituents are `m_i = {1, 1, 1}` and the charges `q_i = {+1, −1, −1}`. We can estimate the ground state of this Coulombic three-body system using 50 Gaussians

### Psudorandom

```@example 1
using Plots, FewBodyPhysics

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]
masses = [1,1,1]
K = [1/2 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(masses)
K_transformed = J * K * J'
w_transformed = [U' * w_list[i] for i in 1:length(w_list)]
Theortical_value = -0.2620050702328

p, Energy = run_simulation(50,:psudorandom, w_transformed, K_transformed)
plot(p)
```
With a difference in energy
```@example 1
@show Theortical_value - Energy
```
### Quasirandom
And similarly for a quasirandom
```@example 2
using Plots, FewBodyPhysics

w_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]
masses = [1,1,1]
K = [1/2 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(masses)
K_transformed = J * K * J'
w_transformed = [U' * w_list[i] for i in 1:length(w_list)]
Theortical_value = -0.2620050702328

p, Energy = run_simulation(50,:quasirandom, w_transformed, K_transformed)
plot(p)
```
With a difference in energy
```@example 2
@show Theortical_value - Energy
```

### Custom system

One of the strengths of this numerical method is that all the matrix elements are analytical. Consider the following example from [Threshold photoproduction of neutral pions off protons in nuclear model with explicit mesons](https://arxiv.org/pdf/2209.12071.pdf). 

```@example
using Optim, FewBodyPhysics, Plots

b = 3.9
S = 41.5

params = [b, S]
masses = [(m_p+m_n)/2, m_π]

energies, Gaussians, eigenvectors, coordinates, masses = run_simulation_nuclear(5,2,5,masses,params)

rmax = 5 * b
rmin = 0.01 * b
start = log(rmin)
stop = log(rmax)
grid_points = range(rmin,rmax,3000)

Φ = zeros(length(grid_points), length(coordinates))

for i in 1:length(coordinates)
    local ϕ = zeros(length(grid_points))
    ϕ_sum = zeros(length(grid_points))
    rs = coordinates[i]
    c = eigenvectors[i]
    A = [1 / (b^2) for b in rs]
    for j in 2:min(length(c), length(A)) 
        ϕ_sum .+= c[j] .* exp.(-(A[j-1]) .* grid_points.^2)
    end
    ϕ .+= ϕ_sum 
    Φ[:, i] = ϕ ./ c[1]
end

r = range(rmin,rmax, length=3000)

plot(r, Φ[:,1], title="Φ(r)", label="Φ(r)",ylabel="Φ",xlabel="r", linewidth=2) #with phase
Φ_prime = diff(Φ[:,1]) ./ diff(r)
plot!(r[1:end-1], Φ_prime, label="Φ'(r)", linewidth=2)
```