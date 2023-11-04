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

function transform_list(α)
    g_new = [Matrix{Float64}(I, 1, 1) * α[i] for i in 1:length(α)]
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

function S_wave(α, K, w = nothing)
    length = length(α)
    α = transform_list(α)
    kinetic = zeros(length, length)
    overlap = zeros(length, length)
    Coulomb = zeros(length, length)
    
    for i in 1:length
        for j in 1:length
            if j <= i
                A = α[i]
                B = α[j]
                M0, trace, Coul = S_elem(A, B, K, w)
                R = inv(A + B)
                overlap[i, j] = M0
                overlap[j, i] = overlap[i, j]
                kinetic[i, j] = 6 * trace * M0
                kinetic[j, i] = kinetic[i, j]
                Coulomb[i, j] = Coul
                Coulomb[j, i] = Coulomb[i, j]
            end
        end
    end
    return overlap, kinetic, Coulomb
end

function energy_S_wave(bij,K::Matrix{Float64},w)
    α = []
    dim = length(w)
    for i in range(1,length(bij),dim)
        A = trace(A_generate(bij[i:1+dim],w))
        α = push!(A)
    end
    N, T = S_wave(α,K,w)
    H = T+Coulomb
    E = eigen(H,N)
    ground_state_energy = minimum(E)
    return ground_state_energy
end

