using LinearAlgebra

using LinearAlgebra

function jacobi_transform(m_list::Vector{T}) where T
    dim = length(m_list)
    J = Matrix{T}(zeros(dim, dim))
    
    for i in 1:dim
        sum_m = sum(m_list[1:i])
        for j in 1:dim
            if j == i + 1
                J[i, j] = -1
            elseif i + 1 < j
                J[i, j] = 0
            else
                J[i, j] = m_list[j] / sum_m
            end
            if isnan(J[i, j])
                J[i, j] = 1
            end
        end
    end
    
    U = inv(J)
    
    if dim > 1
        U = U[:, 1:end-1]
        J = J[1:end-1, :]
    end
    
    return J, U
end

function corput(n::Int, b=3)
    q = 0.0
    bk = 1.0 / b
    while n > 0
        n, rem = divrem(n, b)
        q += rem * bk
        bk /= b
    end
    return q
end

function halton(n::Int, d)
    x = zeros(Float64, d)
    base = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281]
    @assert length(base) >= d "Error: d exceeds the number of basis elements."
    for i in 1:d
        x[i] = corput(n, base[i])
    end
    return x
end

function w_gen_3()
    w_list = [zeros(3) for _ in 1:3]
    w_list[1][1] = 1
    w_list[1][2] = -1
    w_list[2][1] = 1
    w_list[2][3] = -1
    w_list[3][2] = 1
    w_list[3][3] = -1
    return w_list
end
using LinearAlgebra 

#A_{ij}r_i*r_j
function shift(a::Matrix{T}, b::Matrix{T}) where T
    n = size(a, 2)
    total = zero(T)
    mat = I
    for i in 1:n
        for j in 1:n
            dot_product = dot(a[:, i], b[:, j])
            total += mat[i, j] * dot_product
        end
    end
    return total
end

function transform_list(alphas)
    g_new = [Matrix{Float64}(I, 1, 1) * alphas[i] for i in 1:length(alphas)]
    return g_new
end

function S_elements(A, B, K)
    dim = size(A, 1)
    D = A + B
    R = inv(D)
    M0 = (π^dim / det(D))^(3/2)
    trace = tr(B * K * A * R)
    return M0, trace
end

function A_generate(bij, w_list)
    if isa(bij, Array) || isa(bij, AbstractArray)
        dim = length(w_list)
        mat_list = [w_list[i] * w_list[i]' for i in 1:dim]
        for i in 1:dim
            mat_list[i] .= mat_list[i] / (bij[i]^2)
        end
        A = sum(mat_list)
        return A
    else
        dim = length(w_list)
        mat_list = [w_list[i] * w_list[i]' for i in 1:dim]
        for i in 1:dim
            mat_list[i] .= mat_list[i] / (bij^2)
        end
        A = sum(mat_list)
        return A
    end
end



function S_wave(alphas, K, w = nothing)
    length = length(alphas)
    alphas = transform_list(alphas)
    kinetic = zeros(length, length)
    overlap = zeros(length, length)
    coulomb = zeros(length, length)
    
    for i in 1:length
        for j in 1:length
            if j <= i
                A = alphas[i]
                B = alphas[j]
                M0, trace, coul = S_elem(A, B, K, w)
                R = inv(A + B)
                overlap[i, j] = M0
                overlap[j, i] = overlap[i, j]
                kinetic[i, j] = 6 * trace * M0
                kinetic[j, i] = kinetic[i, j]
                coulomb[i, j] = coul
                coulomb[j, i] = coulomb[i, j]
            end
        end
    end
    
    return overlap, kinetic, coulomb
end

function energyS(bij, K, w)
    alphas = []
    dim = length(w)
    
    for i in 1:dim:length(bij)
        A = A_generate(bij[i:i+dim-1], w)
        push!(alphas, A)
    end
    
    N, kinetic, coulomb = S_wave(alphas, K, w)
    H = kinetic + coulomb
    E = eigen(H).values
    E0 = minimum(E)
    
    return E0
end
