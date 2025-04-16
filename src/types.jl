module Types

export Particle,
       GaussianBase, Rank0Gaussian, Rank1Gaussian, Rank2Gaussian,
       BasisSet,
       Operator, KineticEnergy, CoulombPotential,
       FewBodyHamiltonian, MatrixElementResult

abstract type GaussianBase end

struct Particle
    mass::Float64
    charge::Float64
    label::Symbol
end

struct Rank0Gaussian <: GaussianBase
    A::Matrix{Float64}  # Correlation matrix
end

struct Rank1Gaussian <: GaussianBase
    A::Matrix{Float64}                 # Correlation matrix
    a::Vector{Vector{Float64}}         # Polarization vectors
end

struct Rank2Gaussian <: GaussianBase
    A::Matrix{Float64}
    a::Vector{Vector{Float64}}         # First polarization vector set
    b::Vector{Vector{Float64}}         # Second polarization vector set
end

struct BasisSet
    functions::Vector{GaussianBase}
end

abstract type Operator end

struct KineticEnergy <: Operator
    K::Matrix{Float64}  # required for ⟨T⟩
end

struct CoulombPotential <: Operator
    coefficient::Float64
    w::Vector{Float64}
end


struct FewBodyHamiltonian
    basis::BasisSet
    operators::Vector{Operator}
end

struct MatrixElementResult
    bra::GaussianBase
    ket::GaussianBase
    operator::Operator
    value::Float64
end

end # module