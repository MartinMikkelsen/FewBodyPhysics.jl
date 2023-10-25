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