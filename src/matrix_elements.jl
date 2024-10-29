using LinearAlgebra
using Optim

export S_elements, S_wave, S_energy, P_elements, pion_nucleon,ComputeEigenSystem, GetMinimumEnergy, OptimizeGlobalParameters

struct PositiveDefiniteSymmetricMatrix{T<:Real}
    matrix::Matrix{T}
end

function PositiveDefiniteSymmetricMatrix(mat::Matrix{T}) where T <: Real
    @assert size(mat, 1) == size(mat, 2) "Matrix must be square."
    @assert issymmetric(mat) "Matrix must be symmetric."
    @assert isposdef(mat) "Matrix must be positive definite."
    return PositiveDefiniteSymmetricMatrix{T}(mat)
end

"""
    S_elements(A::Matrix{Float64}, B::Matrix{Float64}, K::Matrix{Float64}, w::Union{Nothing, Vector{Vector{Float64}}}=nothing)

Calculate matrix elements for overlap, kinetic energy, and optionally the Coulomb term.

# Arguments
- `A`: Width matrix of Gaussian basis functions for state `i`.
- `B`: Width matrix of Gaussian basis functions for state `j`.
- `K`: Kinetic energy matrix.
- `w` (optional): List of weight vectors for the particles involved.

# Returns
- `M0`: The overlap matrix element between the two states.
- `trace`: The trace used in the kinetic energy calculation.
- `Coulomb_term` (optional): The Coulomb interaction term, if weight vectors `w` are provided.
"""
function S_elements(A::Matrix{Float64}, B::Matrix{Float64}, K::Matrix{Float64}, w::Union{Nothing, Vector{Vector{Float64}}}=nothing)
    dim = size(A, 1)
    D = A + B
    @assert isposdef(D) "Matrix D must be positive definite."
    R = inv(D)
    detD = det(D)
    @assert detD > 0 "Determinant of D must be positive."
    M0 = (π^dim / detD)^(1.5)
    trace = tr(B * K * A * R)
    Coulomb_term = 0.0
    if w !== nothing
        for k in 1:length(w)
            wk = w[k]
            β = 1.0 / (wk' * R * wk)
            factor = ifelse(k == 3, 2.0, -2.0)
            Coulomb_term += factor * sqrt(β / π) * M0
        end
        return M0, trace, Coulomb_term
    else
        return M0, trace
    end
end

"""
    S_wave(α::Vector{Matrix{Float64}}, K::Matrix{Float64}, w::Union{Nothing, Vector{Vector{Float64}}}=nothing)

Calculate the overlap, kinetic energy, and optionally Coulomb interaction matrices for a given set of basis functions.

# Arguments
- `α`: A list of width matrices for the Gaussian basis functions.
- `K`: Kinetic energy matrix.
- `w` (optional): List of weight vectors for the particles involved.

# Returns
- `overlap`: The overlap matrix for the basis functions.
- `kinetic`: The kinetic energy matrix for the basis functions.
- `Coulomb`: The Coulomb interaction matrix, if weight vectors `w` are specified.
"""
function S_wave(α::Vector{Matrix{Float64}}, K::Matrix{Float64}, w::Union{Nothing, Vector{Vector{Float64}}}=nothing)
    len = length(α)
    overlap = zeros(Float64, len, len)
    kinetic = zeros(Float64, len, len)
    Coulomb = w !== nothing ? zeros(Float64, len, len) : nothing
    for i in 1:len
        for j in 1:i
            A, B = α[i], α[j]
            if w !== nothing
                M0, trace, Coulomb_term = S_elements(A, B, K, w)
                Coulomb[i, j] = Coulomb[j, i] = Coulomb_term
            else
                M0, trace = S_elements(A, B, K)
            end
            overlap[i, j] = overlap[j, i] = M0
            kinetic[i, j] = kinetic[j, i] = 6 * trace * M0
        end
    end
    return overlap, kinetic, Coulomb
end

"""
    S_energy(bij::Vector{Float64}, K::Matrix{Float64}, w::Vector{Vector{Float64}})

Compute the ground state energy of the system using the basis functions specified by the width parameters `bij`.

# Arguments
- `bij`: A vector of width parameters for the Gaussian basis functions.
- `K`: Kinetic energy matrix.
- `w`: List of weight vectors for the particles involved.

# Returns
- `E0`: The lowest eigenvalue computed from the Hamiltonian.
"""
function S_energy(bij::Vector{Float64}, K::Matrix{Float64}, w::Vector{Vector{Float64}})
    α = Vector{Matrix{Float64}}()
    dim = length(w)
    @assert length(bij) % dim == 0 "Length of bij must be a multiple of the dimension."
    for i in 1:dim:length(bij)
        bij_segment = bij[i:i+dim-1]
        A = generate_A_matrix(bij_segment, w)
        push!(α, A)
    end
    overlap, kinetic, Coulomb = S_wave(α, K, w)
    H = kinetic + Coulomb
    E, _ = eigen(H, overlap)  
    E0 = minimum(real(E))  
    return E0
end


"""
    P_elements(a::Vector{Float64}, b::Vector{Float64}, A::PositiveDefiniteSymmetricMatrix, B::PositiveDefiniteSymmetricMatrix, K::Matrix{Float64}, w::Union{Nothing, Vector{Vector{Float64}}}=nothing)

Calculate the perturbation matrix elements given two basis states represented by vectors `a` and `b`, and their respective width matrices `A` and `B`.

# Arguments
- `a`: The coefficient vector for basis state `i`.
- `b`: The coefficient vector for basis state `j`.
- `A`: Width matrix for state `i`.
- `B`: Width matrix for state `j`.
- `K`: Kinetic energy matrix.
- `w` (optional): List of weight vectors for the particles involved.

# Returns
- `M1`: The overlap perturbation term.
- `kinetic`: The kinetic energy perturbation term.
- `Coulomb_term` (optional): The Coulomb interaction perturbation term, if weight vectors `w` are provided.
"""
function P_elements(a::Vector{Float64}, b::Vector{Float64}, A::PositiveDefiniteSymmetricMatrix, B::PositiveDefiniteSymmetricMatrix, K::Matrix{Float64}, w::Union{Nothing, Vector{Vector{Float64}}}=nothing)
    D = A.matrix + B.matrix
    @assert isposdef(D) "Matrix D must be positive definite."
    R = inv(D)
    M0, trace = S_elements(A.matrix, B.matrix, K)
    M1 = 0.5 * (a' * R * b) * M0 

    kinetic = 6 * trace * M1
    kinetic += (a' * K * b) * M0
    kinetic -= (a' * K * A.matrix * R * b) * M0
    kinetic -= (a' * R * B.matrix * K * b) * M0
    kinetic += (a' * R * B.matrix * K * A.matrix * R * b) * M0

    if w !== nothing
        w_concat = hcat(w...)
        β = 1.0 / (w_concat' * R * w_concat)[1]
        Coulomb_term = 2 * sqrt(β / π) * M1
        Coulomb_term -= sqrt(β / π) * β / 3 * (a' * R * w_concat * w_concat' * R * b) * M0
        return M1, kinetic, Coulomb_term
    else
        return M1, kinetic
    end
end

"""
    pion_nucleon(alphas, masses, params)

Calculate the overlap and kinetic matrices for a pion-nucleon system.

# Arguments
- `alphas`: A vector of alpha values, which are parameters related to the Gaussian basis functions.
- `masses`: A 2-element vector containing the masses of the nucleon and the pion, respectively.
- `params`: A 2-element vector containing the parameters `b` and `S`.

# Returns
- `overlap`: A matrix representing the overlap between the basis functions.
- `kinetic`: A matrix representing the kinetic energy elements.

# Description
The function calculates the overlap and kinetic matrices for a pion-nucleon system using Gaussian basis functions. The overlap matrix elements are calculated as integrals of the product of two basis functions, and the kinetic matrix elements are calculated as integrals of the product of the derivatives of two basis functions.

The function uses the reduced mass of the pion-nucleon system, the parameter `b` related to the width of the Gaussian functions, and the parameter `S` related to the amplitude of the Gaussian functions.

The matrices are symmetric, and the diagonal elements of the overlap matrix are 1. The off-diagonal elements are calculated using the alpha parameters and the `b` and `S` parameters.
"""
function pion_nucleon(alphas, masses, params)
    ħc = 197.3
    b = params[1]
    S = params[2]
    mass_nucleon = masses[1]
    mass_π = masses[2]
    mass_reduced = mass_nucleon * mass_π / (mass_nucleon + mass_π)
    K = (ħc)^2 / (2 * mass_reduced)
    κ = 1 / b^2
    κ = [κ] 
    Amplitude = S / b 
    length = size(alphas, 1) + 1
    kinetic = zeros(length, length)
    overlap = zeros(length, length) 
    overlap[1, 1] = 1
    kinetic[1, 1] = 0
    for i in 1:length
        for j in 1:i
            if i == 1 && j == 1
                continue
            elseif j == 1 && i != 1
                B = alphas[i - 1]
                overlap[i, j] = 0
                overlap[j, i] = overlap[i, j]
                kinetic[i, j] = 3 * Amplitude * 3 / 2 * 1 / (B + κ[1]) * (π / (B + κ[1]))^(3 / 2)
                kinetic[j, i] = kinetic[i, j]
            else 
                A = alphas[i - 1]; B = alphas[j - 1]
                overlap[i, j] = 3 * 3 / 2 * 1 / (B + A) * (π / (B + A))^(3 / 2)
                overlap[j, i] = overlap[i, j]
                kinetic[i, j] = 3 * K * 15 * A * B / ((A + B)^2) * (π / (A + B))^(3 / 2) + mass_π * overlap[i, j]
                kinetic[j, i] = kinetic[i, j]
            end
        end
    end
    return overlap, kinetic
end

"""
    ComputeEigenSystem(bs, masses, params)

Calculate the eigenvalues and eigenvectors of a system defined by parameters `bs`, `masses`, and `params`.

# Arguments
- `bs`: Array of parameter values used in the computation.
- `masses`: Array of masses, representing physical properties of the system.
- `params`: Additional parameters required for the calculation.

# Returns
- Tuple of eigenvalues (`E`) and eigenvectors (`c`).
"""
function ComputeEigenSystem(bs, masses, params)
    A = [1 / b^2 for b in bs]  
    N, H = pion_nucleon(A, masses, params)
    E, c = eigen(H, N)
    return E, c
end
"""
    GetMinimumEnergy(bs, masses, params)

Compute the minimum energy of a system characterized by `bs`, `masses`, and `params`.

# Arguments
- `bs`: Array of parameter values used in the computation.
- `masses`: Array of masses, representing physical properties of the system.
- `params`: Additional parameters required for the calculation.

# Returns
- Minimum energy value of the system.
"""
function GetMinimumEnergy(bs, masses, params)
    E, _ = ComputeEigenSystem(bs, masses, params)
    return E[1] 
end
"""
    OptimizeGlobalParameters(ngauss, dim, bmax, masses, params)

Perform global optimization over a given parameter space to find optimal parameters for a physical system.

# Arguments
- `ngauss`: Number of Gaussian functions used in the optimization.
- `dim`: Dimensionality of the parameter space.
- `bmax`: Maximum value of the parameter `b`.
- `masses`: Array of masses, representing physical properties of the system.
- `params`: Additional parameters used in the optimization.

# Returns
- `E_list`: List of optimized energies.
- `gaussians`: List of Gaussian functions used.
- `coords`: Optimized coordinates in the parameter space.
- `eigenvectors`: Eigenvectors corresponding to the optimized coordinates.
- `masses`: Updated masses array.
"""
function OptimizeGlobalParameters(ngauss, dim, bmax, masses, params)
    E_list = []
    gaussians = []
    eigenvectors = []
    coords = []

    global E0S = 0.0
    masses_min = copy(masses)
    masses_min = convert(Vector{Float64}, masses_min)  # Convert to floating-point array
    n_calls = 2

    for i in 1:ngauss
        halt = halton(i, dim)
        bs = log.(halt) * bmax
        for j in 1:n_calls
            masses_min[1] = masses[1] .- E0S
            resS = optimize(bs -> GetMinimumEnergy(bs, masses_min, params), bs, NelderMead(), Optim.Options())
            global optimized_bs = Optim.minimizer(resS)  
            global E0S, C0S = ComputeEigenSystem(optimized_bs, masses_min, params)
            global E0S = E0S[1]  
        end
        push!(E_list, E0S)
        push!(coords, optimized_bs)
        push!(eigenvectors, C0S)
        push!(gaussians, i)
    end
    return E_list, gaussians, eigenvectors, coords, masses
end