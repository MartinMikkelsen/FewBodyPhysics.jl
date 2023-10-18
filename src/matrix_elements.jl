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

function transform_list(alphas)
    g_new = [Matrix{Float64}(I, 1, 1) * alphas[i] for i in 1:length(alphas)]
    return g_new
end

function create_shift(dim, p = nothing, d = nothing)
    if p !== nothing
        a = zeros(3, dim)
        for i in 1:dim
            a[1, i] = 1
        end
        a = svector(a)
        return a
    elseif d !== nothing
        a = zeros(3, dim)
        b = zeros(3, dim)
        for i in 1:dim
            a[1, i] = 1
            b[2, i] = 1
        end
        a = svector(a)
        b = svector(b)
        return a, b
    end
end

function create_pion_shift(dim)
    a = zeros(3, dim)
    b = zeros(Complex{Float64}, 3, dim)
    for i in 1:dim
        a[3, i] = 1
        b[1, i] = 1
        b[2, i] = 1im
    end
    bH = adjoint(b)
    a = svector(a)
    b = svector(b)
    bH = svector(bH)
    return a,b, bH
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
        return M0, trace, coul
    else
        return M0, trace
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

function pion_test_1d(alphas, masses, params)
    b_w = params[1]
    S_w = params[2]
    mN = masses[1]
    mpi = masses[2]
    mNpi = mN * mpi / (mN + mpi)
    K = (197.3)^2 / (2 * mNpi)
    kap = 1 / b_w^2
    kap = I * Matrix(1.0I, 1, 1) * kap
    Amp = S_w / b_w  # Initialization
    length = length(alphas) + 1  # Make dimension one greater due to hardcoding parameter
    kinetic = zeros(Complex{Float64}, length, length)
    overlap = zeros(Complex{Float64}, length, length)  # Initialize matrices
    overlap[1, 1] = 1
    kinetic[1, 1] = 0  # Hardcoding parameters
    for i in 1:length
        for j in 1:length
            if j <= i
                if i == 1 && j == 1
                    continue
                elseif j == 1 && i != 1  # Creation elements
                    B = alphas[i - 1]
                    overlap[i, j] = 0
                    overlap[j, i] = overlap[i, j]
                    kinetic[i, j] = (3 * Amp * 3 / 2) / (B + kap) * (π / (B + kap))^(3 / 2)
                    kinetic[j, i] = kinetic[i, j]
                else  # Kinetic terms
                    A = alphas[i - 1]
                    B = alphas[j - 1]
                    overlap[i, j] = (3 * 3 / 2) / (B + A) * (π / (B + A))^(3 / 2)
                    overlap[j, i] = overlap[i, j]
                    kinetic[i, j] = (3 * K * 15 * A * B / (A + B)^2) * (π / (A + B))^(3 / 2) + mpi * overlap[i, j]
                    kinetic[j, i] = kinetic[i, j]
                end
            end
        end
    end
    return overlap, kinetic
end

