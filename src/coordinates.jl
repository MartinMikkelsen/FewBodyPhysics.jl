using LinearAlgebra

function relative(m)
    N = length(m)
    U_Jacobi = Matrix{Float64}(I, N, N)
    denominator = 1.0
    
    for i in 1:N
        for j in 1:N
            if j < i
                denominator *= m[j]
                U_J[i, j] = m[j] / denominator
            elseif j == i
                U_J[i, j] = -1.0
            else
                U_J[i, j] = 0.0
            end
        end
    end
    
    return U_J
end

# Example usage
m = [1, 2, 3]  # Replace with your specific values
U_J = create_U_J(m)

export relative