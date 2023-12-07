using Optim, LinearAlgebra

export S_elements, S_wave, S_energy, P_elements, pion_nucleon, Energy_pion_photodisintegration

struct PositiveDefiniteSymmetricMatrix{T<:Real}
    matrix::Matrix{T}

    function PositiveDefiniteSymmetricMatrix(mat::Matrix{T}) where T <: Real
        if size(mat, 1) != size(mat, 2)
            throw(ArgumentError("The matrix must be square."))
        end
        if !issymmetric(mat)
            throw(ArgumentError("The matrix must be symmetric."))
        end
        if !isposdef(mat)
            throw(ArgumentError("The matrix must be positive definite."))
        end
        return new{T}(mat)
    end
end
"""
    S_elements(A, B, K, w=nothing)

Calculate matrix elements for overlap, kinetic energy, and optionally the Coulomb term.

# Arguments
- `A::Matrix`: Matrix representing the width of Gaussian basis functions for state `i`.
- `B::Matrix`: Matrix representing the width of Gaussian basis functions for state `j`.
- `K::Matrix`: Kinetic energy matrix.
- `w::Vector` (optional): Weight vectors for the particles involved.

# Returns
- `M0::Float64`: The overlap matrix element between the two states.
- `tra::Float64`: The trace used in the kinetic energy calculation.
- `Coulomb_term::Float64` (optional): The Coulomb interaction term, if weight vectors `w` are provided.

# Notes
- The Coulomb term is calculated only if the weight vectors `w` are specified.
"""
S_elements(A, B, K, w=nothing) = begin
    dim = size(A, 1)
    Coulomb_term = 0.0
    D = A + B
    R = inv(D)
    M0 = (π^dim / det(D))^(3.0 / 2)
    trace = tr(B * K * A * R)
    if !isnothing(w)
        for k in 1:length(w)
            β = 1 / (w[k]' * R * w[k])
            Coulomb_term += (k == 3 ? 2 : -2) * sqrt(β / π) * M0
        end
        return M0, trace, Coulomb_term
    else
        return M0, trace
    end
end
"""
    S_wave(α, K, w=nothing)

Calculate the wavefunction overlap, kinetic energy, and optionally Coulomb interaction matrices for a given set of basis functions.

# Arguments
- `α::Vector`: A list of scalar width parameters for the Gaussian basis functions.
- `K::Matrix`: Kinetic energy matrix.
- `w::Vector` (optional): Weight vectors for the particles involved.

# Returns
- `overlap::Matrix`: The overlap matrix for the basis functions.
- `kinetic::Matrix`: The kinetic energy matrix for the basis functions.
- `Coulomb::Matrix` (optional): The Coulomb interaction matrix, if weight vectors `w` are specified.

# Notes
- The Coulomb matrix is computed only if the weight vectors `w` are specified.
"""
S_wave(α, K, w=nothing) = begin
    len = length(α)
    α = transform_list(α)
    overlap = zeros(len, len)
    kinetic = zeros(len, len)
    Coulomb = zeros(len, len)
    for i in 1:len
        for j in 1:i
            A, B = α[i], α[j]
            M0, trace, Coulomb_term = S_elements(A, B, K, w)
            overlap[i, j] = overlap[j, i] = M0
            kinetic[i, j] = kinetic[j, i] = 6 * trace * M0
            Coulomb[i, j] = Coulomb[j, i] = Coulomb_term
        end
    end
    return overlap, kinetic, Coulomb
end
"""
    S_energy(bij, K, w)

Compute the ground state energy of the system using the basis functions specified by the width parameters `bij`.

# Arguments
- `bij::Vector`: A list of width parameters for the Gaussian basis functions.
- `K::Matrix`: Kinetic energy matrix.
- `w::Vector`: Weight vectors for the particles involved.

# Returns
- `E0::Float64`: The lowest eigenvalue computed from the Hamiltonian, considered as the ground state energy of the system.

# Notes
- This function constructs the Hamiltonian from the overlap, kinetic, and Coulomb matrices and solves for its eigenvalues.
"""
function S_energy(bij, K, w)
    α = []
    dim = length(w)
    for i in 1:dim:length(bij)
        A = A_generate(bij[i:i+dim-1], w)
        push!(α, A)
    end
    N, kinetic, Coulomb = S_wave(α, K, w)
    H = kinetic + Coulomb
    E,v = eigen(H, N)
    E0 = minimum(E)
    return E0
end
"""
    P_elements(a, b, A, B, K, w=nothing)

Calculate the perturbation matrix elements given two basis states represented by vectors `a` and `b`, and their respective width matrices `A` and `B`.

# Arguments
- `a::Vector`: The coefficient vector for basis state `i`.
- `b::Vector`: The coefficient vector for basis state `j`.
- `A::Matrix`: Matrix representing the width of Gaussian basis functions for state `i`.
- `B::Matrix`: Matrix representing the width of Gaussian basis functions for state `j`.
- `K::Matrix`: Kinetic energy matrix.
- `w::Vector` (optional): Weight vectors for the particles involved.

# Returns
- `M1::Float64`: The overlap perturbation term.
- `kinetic::Float64`: The kinetic energy perturbation term.
- `Coulomb_term::Float64` (optional): The Coulomb interaction perturbation term, if weight vectors `w` are provided.

# Notes
- The Coulomb interaction perturbation term is calculated only if the weight vectors `w` are specified.
"""
function P_elements(a, b, A::PositiveDefiniteSymmetricMatrix, B, K, w=nothing)
    D = A + B
    R = inv(D)
    M0, trace = S_elements(A, B, K)
    M1 = 1/2 * (a' * R * b) * M0  # Overlap
    kinetic = 6 * trace * M1  # Kinetic matrix elements
    kinetic += (a' * K * b) * M0
    kinetic += -(a' * (K * A * R) * b) * M0 
    kinetic += -(a' * (R * B * K) * b) * M0 
    kinetic += (a' * (R * B * K * A * R) * b) * M0
    kinetic += (a' * (R * B * K * A * R) * b) * M0

    if w !== nothing
        β = 1 / ((w' * R * w)[1])  # Coulomb terms
        Coulomb_term = 2 * sqrt(β / π) * M1
        Coulomb_term += -sqrt(β / π) * β / 3 * (a' * (R * w * w' * R) * b) * M0 
        return M1, kinetic, Coulomb_term
    else
        return M1, kinetic
    end
end

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
    A = [1 / b^2 for b in bs]  #
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
    coords = []
    eigenvectors = []

    global E0S = 0.0
    masses_min = copy(masses)
    n_calls = 2

    for i in 1:ngauss
        halt = halton(i, dim)
        bs = log.(halt) * bmax
        for j in 1:n_calls
            masses_min[1] = masses[1] .- E0S
            println("Iteration of masses: ", j)
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