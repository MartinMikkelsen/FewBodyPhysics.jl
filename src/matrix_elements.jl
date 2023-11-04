shift_dot(a, b, mat=nothing) = begin
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


w_gen(dim, i, j) = dim == 1 ? [1] : [k == i && 1 || k == j && -1 || 0 for k in 1:dim]

S_elem(A, B, K, w=nothing) = begin
    dim = size(A, 1)
    coul = 0.0
    D = A + B
    R = inv(D)
    M0 = (π^dim / det(D))^(3.0 / 2)
    tra = tr(B * K * A * R)
    if !isnothing(w)
        for k in 1:length(w)
            beta = 1 / (w[k]' * R * w[k])
            coul += (k == 3 ? 2 : -2) * sqrt(beta / π) * M0
        end
        return M0, tra, coul
    else
        return M0, tra
    end
end

S_wave(alphas, K, w=nothing) = begin
    len = length(alphas)
    alphas = transform_list(alphas)
    overlap = zeros(len, len)
    kinetic = zeros(len, len)
    coulomb = zeros(len, len)
    for i in 1:len
        for j in 1:i
            A, B = alphas[i], alphas[j]
            M0, trace, coul = S_elem(A, B, K, w)
            overlap[i, j] = overlap[j, i] = M0
            kinetic[i, j] = kinetic[j, i] = 6 * trace * M0
            coulomb[i, j] = coulomb[j, i] = coul
        end
    end
    return overlap, kinetic, coulomb
end

w_list = w_gen_3()
m_list = [1, 1, 1]
K = [0 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(m_list)
K_trans = J * K * J'
w_trans = [U' * w_list[i] for i in 1:length(w_list)]

function energyS(bij, K, w)
    alphas = []
    dim = length(w)
    for i in 1:dim:length(bij)
        A = A_generate(bij[i:i+dim-1], w)
        push!(alphas, A)
    end
    N, kinetic, coulomb = S_wave(alphas, K, w)
    H = kinetic + coulomb
    E,v = eigen(H, N)
    E0 = minimum(E)
    return E0
end