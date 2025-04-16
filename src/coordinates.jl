module Coordinates

using LinearAlgebra
using ..Types

export ParticleSystem, jacobi_transform, generate_A_matrix, transform_list,
       shift_vectors, generate_weight_vector, transform_coordinates, inverse_transform_coordinates

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

function transform_list(α::Vector{Float64})::Vector{Matrix{Float64}}
    return [Matrix{Float64}([α_i]) for α_i in α]
end

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