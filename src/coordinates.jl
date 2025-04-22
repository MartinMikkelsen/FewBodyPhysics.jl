module Coordinates

using LinearAlgebra
using ..Types

export ParticleSystem, jacobi_transform, generate_A_matrix, transform_list,
       shift_vectors, generate_weight_vector, transform_coordinates, inverse_transform_coordinates, default_b0


"""
    struct ParticleSystem

A structure representing a system of particles with associated masses and coordinate transformations.

# Fields
- `masses::Vector{Float64}`: A vector containing the masses of the particles in the system. Must contain at least two elements.
- `J::Matrix{Float64}`: The Jacobi transformation matrix for the particle system, computed based on the masses.
- `U::Matrix{Float64}`: An auxiliary transformation matrix for the particle system, computed based on the masses.
- `scale::Union{Symbol,Nothing}`: An optional symbol indicating the scale of the system (e.g., `:atomic`, `:molecular`, `:nuclear`). Defaults to `nothing`.

# Constructor
- `ParticleSystem(masses::Vector{Float64}; scale::Union{Symbol,Nothing}=nothing)`: 
  Creates a new `ParticleSystem` instance. The `masses` vector must contain at least two elements. The `scale` parameter is optional and can be used to specify the scale of the system. The Jacobi and auxiliary transformation matrices (`J` and `U`) are computed internally using the `jacobi_transform` function.
"""
struct ParticleSystem
    masses::Vector{Float64}
    J::Matrix{Float64}
    U::Matrix{Float64}
    scale::Union{Symbol,Nothing}  # :atomic, :molecular, :nuclear, etc.

    function ParticleSystem(masses::Vector{Float64}; scale::Union{Symbol,Nothing}=nothing)
        @assert length(masses) ≥ 2 "At least two masses are required for a particle system."
        J, U = jacobi_transform(masses)
        new(masses, J, U, scale)
    end
end

"""
    jacobi_transform(masses::Vector{Float64})::Tuple{Matrix{Float64}, Matrix{Float64}}

Computes the Jacobi transformation matrix `J` and its pseudoinverse `U` for a given vector of masses.

# Arguments
- `masses::Vector{Float64}`: A vector of masses for the system. Must contain at least two masses.

# Returns
- `Tuple{Matrix{Float64}, Matrix{Float64}}`: A tuple containing:
  - `J::Matrix{Float64}`: The Jacobi transformation matrix.
  - `U::Matrix{Float64}`: The pseudoinverse of the Jacobi transformation matrix.

# Details
The Jacobi transformation is used to convert the coordinates of a system of particles into a set of relative coordinates. The transformation matrix `J` is constructed based on the masses of the particles, and its pseudoinverse `U` is computed using the `pinv` function.

# Constraints
- The input vector `masses` must have a length of at least 2. An assertion is raised if this condition is not met.

"""
function jacobi_transform(masses::Vector{Float64})::Tuple{Matrix{Float64}, Matrix{Float64}}
    N = length(masses)
    @assert N ≥ 2 "At least two masses are required for Jacobi transformation."
    J = zeros(Float64, N - 1, N)

    for k in 1:N - 1
        mk = masses[k + 1]
        Mk = sum(masses[1:k])
        μk = sqrt(mk * Mk / (mk + Mk))
        for j in 1:N
            if j ≤ k
                J[k, j] = μk * masses[j] / Mk
            elseif j == k + 1
                J[k, j] = -μk
            else
                J[k, j] = 0.0
            end
        end
    end

    U = pinv(J)
    return J, U
end


"""
    default_b0(scale::Union{Symbol,Nothing}) -> Float64

Returns a default value for the parameter `b0` based on the provided `scale`.

# Arguments
- `scale::Union{Symbol,Nothing}`: A symbol representing the scale type or `nothing`.
    - `:atomic`: Returns `1.0`, corresponding to the Bohr radius in atomic units.
    - `:molecular`: Returns `3.0`, representing a typical molecular bond length.
    - `:nuclear`: Returns `0.03`, approximately 1 femtometer in atomic units.
    - `nothing`: Returns `1.0` as a fallback default.
"""
function default_b0(scale::Union{Symbol,Nothing})
    scale === :atomic     && return 1.0       
    scale === :molecular  && return 3.0       
    scale === :nuclear    && return 0.03     
    scale === nothing     && return 10.0    
    error("Unknown scale: $scale")
end

function generate_A_matrix(bij::Vector{Float64}, w_list::Vector{Vector{Float64}})::Matrix{Float64}
    @assert length(bij) == length(w_list) "Length of `bij` and `w_list` must be equal."
    dim = length(w_list[1])
    A = zeros(Float64, dim, dim)
    for i in 1:length(bij)
        w = w_list[i]
        @assert length(w) == dim "All weight vectors must have the same dimension."
        A += (w * w') / (bij[i]^2)
    end
    return A
end

"""
    transform_list(α::Vector{Float64})::Vector{Matrix{Float64}}

Transforms a vector of `Float64` values into a vector of `Matrix{Float64}` objects. 
Each element of the input vector `α` is wrapped into a 1x1 matrix and returned as 
an element of the resulting vector.

# Arguments
- `α::Vector{Float64}`: A vector of `Float64` values to be transformed.

# Returns
- `Vector{Matrix{Float64}}`: A vector where each element is a 1x1 matrix containing 
  the corresponding value from the input vector `α`.
"""
function transform_list(α::Vector{Float64})::Vector{Matrix{Float64}}
    return [Matrix{Float64}([α_i]) for α_i in α]
end

"""
    shift_vectors(a::Matrix{Float64}, b::Matrix{Float64}, mat::Union{Nothing, Matrix{Float64}}=nothing) -> Float64

Compute a weighted sum of dot products between columns of two matrices `a` and `b`, 
optionally using a weighting matrix `mat`.

# Arguments
- `a::Matrix{Float64}`: A matrix where each column represents a vector.
- `b::Matrix{Float64}`: A matrix where each column represents a vector. Must have the same number of columns as `a`.
- `mat::Union{Nothing, Matrix{Float64}}`: An optional square weighting matrix. If `nothing` is provided, the identity matrix is used.

# Returns
- `Float64`: The computed weighted sum of dot products.

# Constraints
- The number of columns in `a` and `b` must be the same.
- If `mat` is provided, it must be a square matrix with dimensions equal to the number of columns in `a` and `b`.
"""
function shift_vectors(a::Matrix{Float64}, b::Matrix{Float64}, mat::Union{Nothing, Matrix{Float64}}=nothing)::Float64
    n = size(a, 2)
    @assert n == size(b, 2) "Matrices `a` and `b` must have the same number of columns."
    mat = mat === nothing ? I(n) : mat
    @assert size(mat) == (n, n) "Matrix `mat` must be square with size equal to number of vectors."

    sum_val = 0.0
    for i in 1:n
        for j in 1:n
            sum_val += mat[i, j] * dot(view(a, :, i), view(b, :, j))
        end
    end
    return sum_val
end

function generate_weight_vector(dim::Int, i::Int, j::Int)::Vector{Int}
    @assert 1 ≤ i ≤ dim "Index `i` must be between 1 and $dim."
    @assert 1 ≤ j ≤ dim "Index `j` must be between 1 and $dim."
    w = zeros(Int, dim)
    w[i] = 1
    w[j] = -1
    return w
end

function transform_coordinates(J::Matrix{Float64}, r::Vector{Float64})::Vector{Float64}
    @assert size(J, 2) == length(r) "Matrix `J` columns must match length of vector `r`."
    return J * r
end

function inverse_transform_coordinates(U::Matrix{Float64}, x::Vector{Float64})::Vector{Float64}
    @assert size(U, 1) == length(x) "Matrix `U` rows must match length of vector `x`."
    return U * x
end

end 