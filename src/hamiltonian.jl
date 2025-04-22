module Hamiltonian

using LinearAlgebra
using ..Types
using ..MatrixElements

export build_overlap_matrix, build_operator_matrix, build_hamiltonian_matrix, solve_generalized_eigenproblem

struct IdentityOperator <: Operator end

function compute_overlap_element(bra::GaussianBase, ket::GaussianBase)
    A, B = bra.A, ket.A
    R = inv(A + B)
    return (π^length(R) / det(A + B))^(3/2)
end

function build_overlap_matrix(basis::BasisSet)
    n = length(basis.functions)
    S = zeros(n, n)
    for i in 1:n, j in 1:i
        val = compute_overlap_element(basis.functions[i], basis.functions[j])
        S[i, j] = S[j, i] = val
    end
    return S
end

function build_operator_matrix(basis::BasisSet, op::Operator)
    n = length(basis.functions)
    H = zeros(n, n)
    for i in 1:n, j in 1:i
        val = compute_matrix_element(basis.functions[i], basis.functions[j], op)
        H[i, j] = H[j, i] = val
    end
    return H
end

function build_hamiltonian_matrix(basis::BasisSet, operators::Vector{Operator})
    H = zeros(length(basis.functions), length(basis.functions))
    for op in operators
        H .+= build_operator_matrix(basis, op)
    end
    return H
end

function solve_generalized_eigenproblem(H::Matrix{Float64}, S::Matrix{Float64})
    vals, vecs = eigen(H, S)
    return real(vals), real(vecs)
end

end # module