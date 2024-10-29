using Plots

export corput, halton, run_simulation, run_simulation_nuclear

# Generate the nth element of the van der Corput sequence in base b
function corput(n, b=3)
    q, bk = 0.0, 1 / b
    while n > 0
        n, rem = divrem(n, b)
        q += rem * bk
        bk /= b
    end
    return q
end

# Generate the nth d-dimensional point in the Halton sequence
function halton(n, d)
    base = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281]
    @assert d <= length(base) "Dimension exceeds available bases."
    return [corput(n, base[i]) for i in 1:d]
end

# Helper function for calculating energies in quasi-random methods
function calculate_energies(num_gauss, w_transformed, K_transformed, b1, method::Symbol)
    E_list, bases, gaussians = Float64[], [], Int[]
    E_low = Inf
    for i in 1:num_gauss
        hal = halton(i, 15 * length(w_transformed))
        bij = Float64.(-log.(hal) * b1)
        best_energy, best_base = E_low, nothing

        for j in 1:length(w_transformed):length(hal)
            base_segment = bij[j:j+length(w_transformed)-1]
            E0 = S_energy(base_segment, K_transformed, w_transformed)
            if E0 < best_energy
                best_energy = E0
                best_base = copy(base_segment)
            end
        end

        E_low = min(E_low, best_energy)
        push!(E_list, E_low)
        push!(bases, best_base)
        push!(gaussians, i)
        println(E_low)
    end
    return E_list, bases, gaussians, E_low
end

# Main simulation function
function run_simulation(num_gauss::Int, method::Symbol, w_transformed::Vector{Vector{Float64}}, K_transformed::Matrix{Float64}, plot_result::Bool=true)
    b1 = 10.0

    # Choose appropriate method
    if method == :quasirandom
        println("---------QUASI-RANDOM METHOD---------")
        E_list, bases, gaussians, E_low = calculate_energies(num_gauss, w_transformed, K_transformed, b1, method)

    elseif method == :quasirandomrefined
        println("---------QUASI-RANDOM METHOD WITH REFINEMENT---------")
        E_list, bases, gaussians, E_low = calculate_energies(num_gauss, w_transformed, K_transformed, b1, method)

        # Additional refinement step
        for i in 1:length(bases)-length(w_transformed)
            rand_ref = Float64.(-log.(rand(200 * length(w_transformed))) * b1)
            for j in 1:length(w_transformed):length(rand_ref)
                refined_base = bases[i][1:end-length(w_transformed)]  # Adjust size accordingly
                refined_base .= rand_ref[j:j+length(w_transformed)-1]
                E_refined = S_energy(refined_base, K_transformed, w_transformed)
                if E_refined < E_low
                    E_low = E_refined
                    bases[i] .= refined_base
                end
            end
            println("Refined Energy: ", E_low)
        end

    elseif method == :psudorandom
        println("---------PSEUDO-RANDOM METHOD---------")
        E_list, bases, gaussians, E_low = calculate_energies(num_gauss, w_transformed, K_transformed, b1, method)

    else
        error("Invalid method provided!")
    end

    # Plot results
    println("Best convergent numerical value: ", E_list[end])
    p = plot(gaussians, E_list, marker=:circle, label="Numerical result", linewidth=2)
    title!(p, "S-wave Convergence")
    xlabel!(p, "Number of Gaussians")
    ylabel!(p, "Energy [Hartree]")

    if plot_result
        display(p)
    end

    return p, E_list[end], bases
end

"""
Run a nuclear simulation and print the final energy.
"""
function run_simulation_nuclear(ngauss, dim, bmax, masses, params)
    return OptimizeGlobalParameters(ngauss, dim, bmax, masses, params)
end

