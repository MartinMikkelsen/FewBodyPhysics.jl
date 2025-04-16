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

