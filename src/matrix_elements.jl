module MatrixElements

using LinearAlgebra
using ..Types

export compute_matrix_element

"""
compute_matrix_element(bra, ket, op)

Compute the matrix element ⟨bra|op|ket⟩ using analytic expressions.
"""

# ----------- RANK 0 -----------

function compute_matrix_element(bra::Rank0Gaussian, ket::Rank0Gaussian, op::KineticEnergy)
    A, B = bra.A, ket.A
    K = op.K
    R = inv(A + B)
    M0 = (π^length(R) / det(A + B))^(3/2)
    return 6 * tr(B * K * A * R) * M0
end

function compute_matrix_element(bra::Rank0Gaussian, ket::Rank0Gaussian, op::CoulombPotential)
    A, B, w = bra.A, ket.A, op.w
    R = inv(A + B)
    β = 1 / (dot(w, R * w))
    M0 = (π^length(R) / det(A + B))^(3/2)
    return 2 * sqrt(β / π) * M0
end

# ----------- RANK 1 -----------

function compute_matrix_element(bra::Rank1Gaussian, ket::Rank1Gaussian, op::CoulombPotential)
    A, B, a, b, w = bra.A, ket.A, bra.a, ket.a, op.w
    R = inv(A + B)
    β = 1 / (dot(w, R * w))
    M0 = (π^length(R) / det(A + B))^(3/2)
    M1 = 0.5 * dot(b, R * a) * M0
    q2 = 0.25 * dot(a .+ b, R * (w * w') * (a .+ b))
    return 2 * sqrt(β / π) * M1 - sqrt(β^3 / π) / 3 * q2 * M0
end

function compute_matrix_element(bra::Rank1Gaussian, ket::Rank1Gaussian, op::KineticEnergy)
    A, B, a, b = bra.A, ket.A, bra.a, ket.a
    R = inv(A + B)
    M0 = (π^length(R) / det(A + B))^(3/2)
    M1 = 0.5 * dot(b, R * a) * M0
    prefactor = 1.0  # placeholder for ħ² / 2μ

    T1 = 6 * tr(B * A * R) * M1
    T2 = dot(b, a) * M0
    T3 = dot(a, R * B * A * R * b) * M0
    T4 = dot(b, R * B * a) * M0
    T5 = dot(a, R * A * b) * M0

    return prefactor * (T1 + T2 + T3 + T3 - T4 - T5)
end

end # module