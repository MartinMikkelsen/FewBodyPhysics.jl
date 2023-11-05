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
export S_elements, S_wave, S_energy, P_elements
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