using LinearAlgebra

export ParticleSystem, jacobi_transform, generate_A_matrix, transform_list, shift_vectors, generate_weight_vector, transform_coordinates, inverse_transform_coordinates

"""
    ParticleSystem(masses::Vector{Float64})

A data structure representing a system of particles, storing their masses and associated Jacobi transformation matrices for coordinate transformations.

# Fields
- `masses::Vector{Float64}`: A vector containing the masses of the particles.
- `J::Matrix{Float64}`: The Jacobi transformation matrix, used to convert particle coordinates into Jacobi coordinates.
- `U::Matrix{Float64}`: The pseudoinverse of the Jacobi transformation matrix `J`, used to transform Jacobi coordinates back to the particle coordinate system.

# Constructor
- `ParticleSystem(masses::Vector{Float64})`: Constructs a new `ParticleSystem` instance.
  - **Arguments**:
    - `masses`: A vector of particle masses. At least two masses are required.
"""
struct ParticleSystem
    masses::Vector{Float64}
    J::Matrix{Float64}
    U::Matrix{Float64}

    function ParticleSystem(masses::Vector{Float64})
        @assert length(masses) ≥ 2 "At least two masses are required for a particle system."
        J, U = jacobi_transform(masses)
        new(masses, J, U)
    end
end

"""
    jacobi_transform(masses::Vector{Float64})::Tuple{Matrix{Float64}, Matrix{Float64}}

Compute the Jacobi transformation matrices `J` and `U` for a system of particles with specified masses.

# Arguments
- `masses::Vector{Float64}`: A vector of masses for the particles.

# Returns
- `(J::Matrix{Float64}, U::Matrix{Float64})`: The Jacobi transformation matrix and its pseudoinverse.

# Notes
- The matrices `J` and `U` are used to transform between particle coordinates and Jacobi coordinates.
- The pseudoinverse `U` is used instead of the inverse to handle cases where `J` is not square.
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
    generate_A_matrix(bij::Vector{Float64}, w_list::Vector{Vector{Float64}})::Matrix{Float64}

Generate the matrix `A` for Gaussian basis functions given width parameters `bij` and weight vectors `w_list`.

# Arguments
- `bij::Vector{Float64}`: A vector of width parameters for the Gaussian basis functions.
- `w_list::Vector{Vector{Float64}}`: A list of weight vectors.

# Returns
- `A::Matrix{Float64}`: The sum of weighted outer products of `w_list`, scaled by `bij`.

# Notes
- This function constructs the `A` matrix used in the correlated Gaussian method.
"""
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

Transform a list of scalar values `α` into a list of 1x1 matrices.

# Arguments
- `α::Vector{Float64}`: A list of scalar values.

# Returns
- `Array{Matrix{Float64}}`: A list of 1x1 matrices where each matrix contains one of the scalar values from `α`.
"""
function transform_list(α::Vector{Float64})::Vector{Matrix{Float64}}
    return [Matrix{Float64}([α_i]) for α_i in α]
end

"""
    shift_vectors(a::Matrix{Float64}, b::Matrix{Float64}, mat::Union{Nothing, Matrix{Float64}}=nothing)::Float64

Calculate the weighted sum of the element-wise product of vectors `a` and `b` using matrix `mat`.

# Arguments
- `a::Matrix{Float64}`: A matrix where each column is a vector `a_i`.
- `b::Matrix{Float64}`: A matrix where each column is a vector `b_j`.
- `mat::Union{Nothing, Matrix{Float64}}`: An optional matrix to weight the product (default is the identity matrix).

# Returns
- `sum_val::Float64`: The weighted sum of products.

# Notes
- The matrices `a` and `b` should have the same dimensions.
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

"""
    generate_weight_vector(dim::Int, i::Int, j::Int)::Vector{Int}

Generate a weight vector for the `i`-th and `j`-th coordinates in a space of dimension `dim`.

# Arguments
- `dim::Int`: The dimension of the space.
- `i::Int`: The index for the positive element in the weight vector.
- `j::Int`: The index for the negative element in the weight vector.

# Returns
- `w::Vector{Int}`: A vector with 1 at the `i`-th position, -1 at the `j`-th position, and 0 elsewhere.
"""
function generate_weight_vector(dim::Int, i::Int, j::Int)::Vector{Int}
    @assert 1 ≤ i ≤ dim "Index `i` must be between 1 and $dim."
    @assert 1 ≤ j ≤ dim "Index `j` must be between 1 and $dim."
    w = zeros(Int, dim)
    w[i] = 1
    w[j] = -1
    return w
end

"""
    transform_coordinates(J::Matrix{Float64}, r::Vector{Float64})::Vector{Float64}

Transform the coordinates `r` of a system using the Jacobi matrix `J`.

# Arguments
- `J::Matrix{Float64}`: The Jacobi transformation matrix.
- `r::Vector{Float64}`: The coordinates to be transformed.

# Returns
- `x::Vector{Float64}`: The transformed coordinates.
"""
function transform_coordinates(J::Matrix{Float64}, r::Vector{Float64})::Vector{Float64}
    @assert size(J, 2) == length(r) "Matrix `J` columns must match length of vector `r`."
    return J * r
end

"""
    inverse_transform_coordinates(U::Matrix{Float64}, x::Vector{Float64})::Vector{Float64}

Transform the coordinates `x` back to the original system using the inverse of the Jacobi matrix.

# Arguments
- `U::Matrix{Float64}`: The inverse Jacobi transformation matrix.
- `x::Vector{Float64}`: The coordinates in Jacobi space.

# Returns
- `r::Vector{Float64}`: The coordinates transformed back to the original system.
"""
function inverse_transform_coordinates(U::Matrix{Float64}, x::Vector{Float64})::Vector{Float64}
    @assert size(U, 1) == length(x) "Matrix `U` rows must match length of vector `x`."
    return U * x
end
