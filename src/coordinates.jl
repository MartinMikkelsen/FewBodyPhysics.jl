using LinearAlgebra

export Ω, A_generate, transform_list, shift, w_gen, transform_coordinates, transform_back
"""
    Ω(masses::Array)

Calculate the Jacobi transformation matrix `J` and its inverse `U` for a system of particles with specified `masses`.

# Arguments
- `masses::Array`: A vector of masses for the particles.

# Returns
- `J::Matrix`: The Jacobi transformation matrix.
- `U::Matrix`: The inverse of the Jacobi transformation matrix.

# Notes
- For systems with more than one particle, the returned matrices exclude the last row/column for proper dimensionality in transformations.
"""
Ω(masses::Array) = begin
    dim = size(masses, 1)
    J = zeros(dim, dim)
    for i in 1:dim
        sum_m = sum(masses[1:i])
        for j in 1:dim
            J[i, j] = j == i + 1 ? -1 : i + 1 < j ? 0 : masses[j] / sum_m
            J[i, j] = isnan(J[i, j]) ? 1 : J[i, j]
        end
    end
    U = inv(J)
    U = dim > 1 ? U[:, 1:end-1] : U
    J = dim > 1 ? J[1:end-1, :] : J
    return J, U
end
"""
    A_generate(bij::Array, w_list::Array)

Generate a matrix A for Gaussian basis functions given width parameters `bij` and weight vectors `w_list`.

# Arguments
- `bij::Array`: A vector of width parameters for the Gaussian basis functions.
- `w_list::Array`: A list of weight vectors.

# Returns
- `Matrix`: The sum of weighted outer products of `w_list`, scaled by `bij`.

# Notes
- This function is used to construct basis elements for the expansion of few-body wavefunctions.
"""
A_generate(bij, w_list::Array) = begin
    dim = length(w_list)
    mat_list = [(w_list[i] * w_list[i]') ./ (bij[i]^2) for i in 1:dim]
    return sum(mat_list)
end

"""
    transform_list(α::Array)

Transform a list of scalar values `α` into a list of 1x1 matrices.

# Arguments
- `α::Array`: A list of scalar values.

# Returns
- `Array`: A list of 1x1 matrices where each matrix contains one of the scalar values from `α`.
"""
function transform_list(α::Array)
    return[ones(1, 1) .* α[i] for i in 1:length(α)]
end
"""
    shift(a::Array, b::Array, mat::Matrix=I)

Calculate the weighted sum of the element-wise product of vectors `a` and `b` using matrix `mat`.

# Arguments
- `a::Array`: A vector or matrix.
- `b::Array`: A vector or matrix of the same size as `a`.
- `mat::Matrix`: An optional matrix to weight the product (default is the identity matrix).

# Returns
- `Float64`: The weighted sum of products.

# Notes
- `a` and `b` are typically shift vectors in the configuration space of a few-body system.
- If `mat` is provided, its dimensions must match the number of elements in `a` and `b`.
"""
shift(a, b, mat=nothing) = begin
    n = size(a, 2)
    sum_val = 0.0
    mat = isnothing(mat) ? I : mat
    @assert n == size(mat, 1) "ERROR! Matrix shape does not match number of shift vectors."
    for i in 1:n
        for j in 1:n
            sum_val += mat[i, j] * (a[:, i]' * b[:, j])
        end
    end
    return sum_val
end
"""
    w_gen(dim::Int, i::Int, j::Int)

Generate a weight vector for the i-th and j-th coordinates in a space of dimension `dim`.

# Arguments
- `dim::Int`: The dimension of the space.
- `i::Int`: The index for the positive element in the weight vector.
- `j::Int`: The index for the negative element in the weight vector.

# Returns
- `Vector{Int}`: A vector with 1 at the i-th position, -1 at the j-th position, and 0 elsewhere.

# Notes
- This function is useful for generating basis vectors in few-body coordinate transformations.
"""
w_gen(dim, i, j) = dim == 1 ? [1] : [k == i && 1 || k == j && -1 || 0 for k in 1:dim]

"""
    transform_coordinates(Ω::Matrix{Float64}, r::Vector{Float64})

Transform the coordinates `r` of a system using the Jacobi matrix `Ω`.

# Arguments
- `Ω::Matrix{Float64}`: The Jacobi transformation matrix.
- `r::Vector{Float64}`: The coordinates to be transformed.

# Returns
- `Vector{Float64}`: The transformed coordinates.

# Notes
- This function applies the inverse of Jacobi matrix `J` to the coordinate vector `r`.
"""
function transform_coordinates(Ω::Matrix{Float64}, r::Vector{Float64})
    J, U = Ω(masses)
    return J \ r
end
"""
    transform_back(Ω::Matrix{Float64}, x::Matrix{Float64})

Transform the coordinates `x` back to the original system using the inverse of the Jacobi matrix `Ω`.

# Arguments
- `Ω::Matrix{Float64}`: The Jacobi transformation matrix.
- `x::Matrix{Float64}`: The coordinates to be transformed back.

# Returns
- `Matrix{Float64}`: The coordinates transformed back to the original system.

# Notes
- This function applies the inverse of matrix `U` to the coordinate matrix `x`.
"""
function transform_back(Ω::Matrix{Float64},x::Matrix{Float64})
    J, U = Ω(masses)
    return U \ x
end
