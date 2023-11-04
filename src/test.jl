using LinearAlgebra
using Optim
using Plots
using Random
using PyCall
using Arpack

function jacobi_transform(m_list)
    dim = size(m_list, 1)
    J = zeros(dim, dim)
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

function w_gen_3()
    w1 = zeros(3)
    w2 = zeros(3)
    w3 = zeros(3)
    w1[1] = 1
    w1[2] = -1
    w2[1] = 1
    w2[3] = -1
    w3[2] = 1
    w3[3] = -1
    w_list = [w1, w2, w3]
    return w_list
end

function corput(n, b=3)
    q = 0
    bk = 1 / b
    while n > 0
        n, rem = divrem(n, b)
        q += rem * bk
        bk /= b
    end
    return q
end

function halton(n, d)
    x = zeros(d)
    base = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,233,239,241,251,257,263,269,271,277,281]
    @assert length(base) > d "Error: d exceeds the number of basis elements."
    for i in 1:d
        x[i] = corput(n, base[i])
    end
    return x
end

function A_generate(bij, w_list)
    if isa(bij, AbstractArray)
        dim = length(w_list)
        mat_list = [w_list[i] * w_list[i]' for i in 1:dim]
        for i in 1:dim
            mat_list[i] .= mat_list[i] ./ (bij[i]^2)
        end
        A = sum(mat_list)
        return A
    else
        dim = length(w_list)
        mat_list = [w_list[i] * w_list[i]' for i in 1:dim]
        for i in 1:dim
            mat_list[i] .= mat_list[i] ./ (bij^2)
        end
        A = sum(mat_list)
        return A
    end
end

function shift_dot(a, b, mat=nothing)
    n = size(a, 2)
    sum_val = 0
    if mat === nothing
        mat = I
    end
    @assert n == size(mat, 1) "ERROR! Matrix shape does not match number of shift vectors."
    for i in 1:n
        for j in 1:n
            dot = a[:, i]' * b[:, j]
            sum_val += mat[i, j] * dot
        end
    end
    return sum_val
end

function transform_list(alphas)
    g_new = [ones(1, 1) .* alphas[i] for i in 1:length(alphas)]
    println("size=", size(g_new))
    return g_new
end

function w_gen(dim, i, j)
    if dim == 1
        w = ones(1, 1)
        return w
    else
        w = zeros(1, dim)
        w[i] = 1
        w[j] = -1
        return w
    end
end

np = pyimport("numpy")
sci = pyimport("scipy")
function S_elem(A, B, K, w=nothing)
    dim = size(A, 1)
    coul = 0.0
    D = A + B
    R = inv(D)
    M0 = (π^dim / det(D))^(3.0 / 2)
    trace = np.trace(B * K * A * R)
    if w !== nothing
        for k in 1:length(w)
            beta = 1 / (w[k]' * R * w[k])
            if k == 3
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

function transform_list(alphas)
    g_new = [ones(1, 1) .* alphas[i] for i in 1:length(alphas)]
    return g_new
end

function S_wave(alphas, K, w=nothing)
    len = length(alphas)
    alphas = transform_list(alphas)
    kinetic = zeros(len, len)
    overlap = zeros(len, len)
    coulomb = zeros(len, len)
    for i in 1:len
        for j in 1:len
            if j <= i
                A = alphas[i]
                B = alphas[j]
                M0, trace, coul = S_elem(A, B, K, w)
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

w_list = w_gen_3()
m_list = np.array([1, 1, 1])
K = np.array([[0,0,0],[0,1/2,0],[0,0,1/2]])
J, U = jacobi_transform(m_list)
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
    E = sci.linalg.eigh(H, N, eigvals_only=true)
    E0 = minimum(E)
    return E0
end

function run_simulation()
    b1 = 7
    E_list = []
    gaussians = []
    E_theoS = []
    bij = []
    E_S = -0.527

    println("---------QUASI-RANDOM METHOD---------")

    E_low = Inf
    bases = []
    base_test = []
    for i in 1:50
        hal = halton(i, 15 * length(w_trans))
        bij = -log.(hal) .* b1
        for j in 1:length(w_trans):length(hal)
            append!(base_test, bij[j:j+length(w_trans)-1])
            E0 = energyS(base_test, K_trans, w_trans)
            if E0 <= E_low
                E_low = E0
                global base_curr = copy(bij[j:j+length(w_trans)-1])
            end
            base_test = base_test[1:end-length(w_trans)]
        end
        append!(bases, base_curr)
        append!(base_test, base_curr)
        push!(E_list, E_low)
        println(E_low)
        push!(gaussians, i)
        push!(E_theoS, E_S)
    end

    println("Best convergent numerical value: ", E_list[end])
    println("Theoretical value: ", E_S)
    println("Difference: ", abs(E_list[end] - E_S))

    plot(gaussians, E_list, marker=:circle, label="Numerical result")
    plot!(gaussians, E_theoS, linestyle=:dash, label="Theoretical value")
    title!("S-wave convergence of Positron and two Electron System")
    xlabel!("Number of Gaussians")
    ylabel!("Energy [Hartree]")
end

run_simulation()