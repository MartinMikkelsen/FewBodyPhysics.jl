# FewBodyPhysics

[![Build Status](https://github.com/MartinMikkelsen/FewBodyPhysics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MartinMikkelsen/FewBodyPhysics.jl/actions/workflows/CI.yml?query=branch%3Amain)


Welcome to FewBodyPhysics.jl, a Julia package dedicated to the study and simulation of quantum mechanical few-body systems using explicitly correlated Gaussian methods. This package offers a powerful computational framework to model and analyze various quantum systems, from atoms and molecules to light nuclei and quantum dots.


## Features

- Variational Method Implementation: Utilizes the variational principle in quantum mechanics for approximating ground state energies and other properties.
Correlated Gaussian Basis Functions: Employs explicitly correlated Gaussian functions to accurately represent particle correlations.
- Wide Applicability: Suitable for a range of systems, including small atoms, molecules, and exotic quantum states.
- High Precision: Offers detailed and precise modeling capabilities.

## Basic Theortical Foundations

### Variational Principle
The package is based on the variational principle in quantum mechanics, which can be expressed as:

$$E_0 \leq \frac{\langle \Psi | H | \Psi \rangle}{\langle \Psi | \Psi \rangle},$$
where `E_0` represents the ground state energy, `\Psi` is the trail wave function and `H` is the Hamiltonian of the system.

### Explicitly Correlated Gaussians (ECG)

The ECG function for an N-particle system is given by:

$$\psi_G = \exp\big(-\frac{1}{2} \sum_{i<j}^{N} \alpha_{ij} (r_i - r_j)^2\big),$$
where `r_i` and `r_j` are the position vectors of particles `i` and `j`, respectively, while `\alpha_{ij}` are the variational parameters. 
To use this package, import it in your Julia script:


## Install

To install FewBodyPhysics.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("FewBodyPhysics")
```