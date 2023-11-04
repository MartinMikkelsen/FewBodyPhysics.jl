using LinearAlgebra

Ω(m_list) = begin
    dim = size(m_list, 1)
    J = zeros(dim, dim)
    for i in 1:dim
        sum_m = sum(m_list[1:i])
        for j in 1:dim
            J[i, j] = j == i + 1 ? -1 : i + 1 < j ? 0 : m_list[j] / sum_m
            J[i, j] = isnan(J[i, j]) ? 1 : J[i, j]
        end
    end
    U = inv(J)
    U = dim > 1 ? U[:, 1:end-1] : U
    J = dim > 1 ? J[1:end-1, :] : J
    return J, U
end

w_gen_3() = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]

corput(n, b=3) = begin
    q, bk = 0.0, 1 / b
    while n > 0
        n, rem = divrem(n, b)
        q += rem * bk
        bk /= b
    end
    return q
end

halton(n, d) = begin
    base = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,233,239,241,251,257,263,269,271,277,281]
    @assert length(base) ≥ d "Error: d exceeds the number of basis elements."
    return [corput(n, base[i]) for i in 1:d]
end

A_generate(bij, w_list) = begin
    dim = length(w_list)
    mat_list = [(w_list[i] * w_list[i]') ./ (bij[i]^2) for i in 1:dim]
    return sum(mat_list)
end

function transform_list(alphas)
    return[ones(1, 1) .* alphas[i] for i in 1:length(alphas)]
end

function Λ(Ω::Matrix{T}, m_list::Vector{T}) where T
    dim = length(m_list)
    Λ = Matrix{T}(zeros(dim, dim))
    J, U = Ω(m_list)
    for k in 1:dim
        for j in 1:dim-1
            for i in 1:dim-1
                Λ[i, k] += J[i, k] * J[j, k] / m_list[k]
            end
        end
    end
    return Λ
end
### Consider algorithm
function transform_coordinates(Ω::Matrix{Float64}, r::Vector{Float64})
    J, U = Ω(m_list)
    return J \ r
end

function transform_back(Ω::Matrix{Float64},x::Matrix{Float64})
    J, U = Ω(m_list)
    return U \ x
end
