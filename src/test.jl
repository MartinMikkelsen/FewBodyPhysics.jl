using LinearAlgebra

function Ω(m_list::Vector{T}) where T
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

function energy_S_wave(bij, K::Matrix{Float64}, w)
    α = Vector{Float64}()
    dim = length(w)
    for i in 1:dim:length(bij)
        A = A_generate(bij[1:i], w)
        push!(α, tr(A))  # Corrected trace calculation
    end
    N, T, _ = S_wave(α, K, w)  # I've removed Coulomb term as it was not defined
    H = T
    E = eigen(H, N)
    ground_state_energy = minimum(E.values)
    return ground_state_energy
end

w_list = w_gen_3()
m_list = [Inf, 1.0, 1.0]
K = [0 0 0; 0 1/2 0; 0 0 1/2]
J, U = Ω(m_list)
K_trans = J * K * J'
w_trans = [U' * w_list[i] for i in 1:length(w_list)]

b1=7
E_list=[]
gaussians=[]
E_theoS=[]
bij = Float64[]

E_S=-0.527

E_low = Inf
bases = []
base_test = []
E_list = []
gaussians = Int[]
E_theoS = []

for i in 1:25
    hal = halton(i, 15 * length(w_trans))
    bij = -log.(hal) * b1
    base_curr = Float64[]

    for j in 1:length(hal):length(w_trans)
        push!(base_test, bij[j:j+length(w_trans)-1]...)  # Modify the content of base_test

        E0 = energy_S_wave(base_test, K_trans, w_trans)

        if E0 <= E_low
            E_low = E0
            base_curr = copy(bij[j:j+length(w_trans)-1])
        end

        # Remove elements from base_test
        base_test = base_test[1:end-length(w_trans)]
    end

    append!(bases, base_curr)
    append!(base_test, base_curr)
    push!(E_list, E_low)
    push!(gaussians, i)
    push!(E_theoS, E_S)

    println("E_low for iteration $i: $E_low")
end

println("Best convergent numerical value: ", E_list[end])
println("Theoretical value: ", E_S)
println("Difference: ", abs(E_list[end] - E_S))