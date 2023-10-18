using LinearAlgebra 

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

function S_elements(A::AbstractMatrix{T}, B::AbstractMatrix{T}, K::AbstractMatrix{T}, w::AbstractVector{T}=nothing) where T
    dim = size(A, 1)
    coul = 0.0
    D = A + B
    R = inv(D)
    @show det(D)
    M0 = (π^dim / det(D))^(3/2)
    trace = tr(B * K * A * R)
    
    if w !== nothing
        for k in 1:length(w)
            beta = 1 / (w[k]' * R * w[k])
            if k == 2
                coul += 2 * sqrt(beta / π) * M0
            else
                coul -= 2 * sqrt(beta / π) * M0
            end
        end
        return (M0, trace, coul)
    else
        return (M0, trace)
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
