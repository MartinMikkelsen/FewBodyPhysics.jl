using Plots

export corput, halton, run_simulation, run_simulation_nuclear

"""
Generate the nth element of the van der Corput sequence in base b.

# Arguments
- `n`: The nth element of the sequence to be generated.
- `b`: The base for the van der Corput sequence. Default is 3.

# Returns
- `q`: The nth element of the van der Corput sequence in base b.
"""
corput(n, b=3) = begin
    q, bk = 0.0, 1 / b
    while n > 0
        n, rem = divrem(n, b)
        q += rem * bk
        bk /= b
    end
    return q
end
"""
Generate the nth d-dimensional point in the Halton sequence.

# Arguments
- `n`: The nth element of the sequence to be generated.
- `d`: The dimension of the space.

# Returns
- An array containing the nth d-dimensional point in the Halton sequence.

# Errors
- Throws an assertion error if `d` exceeds the number of basis elements.
"""
halton(n, d) = begin
    base = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,233,239,241,251,257,263,269,271,277,281]
    @assert length(base) ≥ d "Error: d exceeds the number of basis elements."
    return [corput(n, base[i]) for i in 1:d]
end

"""
Run the simulation for a quantum system using quasi-random or pseudo-random methods to determine the S-wave convergence.

# Arguments
- `num_gauss::Int`: The number of Gaussians to use in the simulation. Default is 15.
- `method::Symbol`: The method to use for the simulation. Can be `:quasirandom`, `:quasirandomrefined`, or `:psudorandom`. Default is `:quasirandom`.
- `plot_result::Bool`: Whether to plot the results. Default is true.

# Returns
- `p`: The plot object if `plot_result` is true.

# Notes
- The function prints various convergence information and, if `plot_result` is true, displays a plot of the numerical result against the theoretical value.
"""
function run_simulation(num_gauss::Int, method::Symbol, w_transformed::Array, K_transformed::Array, plot_result::Bool=true)
    b1 = 10
    E_list = []
    gaussians = []
    bij = []
    E_low = 1e10
    global bases = []
    global base_test = []
    if method == :quasirandom
        println("---------QUASI-RANDOM METHOD---------")

        for i in 1:num_gauss
            hal = halton(i, 15 * length(w_transformed))
            bij = -log.(hal) .* b1
            for j in 1:length(w_transformed):length(hal)
                append!(base_test, bij[j:j+length(w_transformed)-1])
                E0 = S_energy(base_test, K_transformed, w_transformed)
                if E0 <= E_low
                    E_low = E0
                    global base_curr = copy(bij[j:j+length(w_transformed)-1])
                end
                base_test = base_test[1:end-length(w_transformed)]
            end
            append!(bases, base_curr)
            append!(base_test, base_curr)
            push!(E_list, E_low)
            println(E_low)
            push!(gaussians, i)
        end

    elseif method == :quasirandomrefined
        println("---------QUASI-RANDOM METHOD WITH REFINEMENT---------")

        for i in 1:num_gauss
            hal = halton(i, 15 * length(w_transformed))
            bij = -log.(hal) .* b1
            for j in 1:length(w_transformed):length(hal)
                append!(base_test, bij[j:j+length(w_transformed)-1])
                E0 = S_energy(base_test, K_transformed, w_transformed)
                if E0 <= E_low
                    E_low = E0
                    global base_curr = copy(bij[j:j+length(w_transformed)-1])
                end
                base_test = base_test[1:end-length(w_transformed)]
            end
            append!(bases, base_curr)
            append!(base_test, base_curr)
            push!(E_list, E_low)
            println(E_low)
            push!(gaussians, i)
        end

        bases_ref = copy(bases)
        E_ref = E_list[end]
        E_list_ref = []
        for i in 1:length(bases_ref)-length(w_transformed)
            rand_ref = rand(200 * length(w_transformed))
            bij_ref = .-log.(rand_ref) * b1
            for j in 1:length(w_transformed):length(rand_ref)
                bases_ref[i:i+length(w_transformed)-1] = bij_ref[j:j+length(w_transformed)-1]
                E_test = S_energy(bases_ref, K_transformed, w_transformed)
                if E_test < E_ref
                    E_ref = E_test
                    bases[i:i+length(w_transformed)-1] = bij_ref[j:j+length(w_transformed)-1]
                end
            end
            bases_ref = copy(bases)
            push!(E_list_ref, E_ref)
            println("E_ref: ", E_ref)
        end

        println("Energy after refinement: ", E_ref)
        println("Difference in energy from before refinement: ", abs(E_ref - E_list[end]))

    elseif method == :psudorandom
        println("---------PSEUDO-RANDOM METHOD---------")
        for i in 1:num_gauss
            rnd = rand(400 * length(w_transformed))
            bij2 = -log.(rnd) * b1
            for j in 1:length(w_transformed):length(rnd)
                append!(base_test, bij2[j:j+length(w_transformed)-1])
                E0 = S_energy(base_test, K_transformed, w_transformed)
                if E0 <= E_low
                    E_low = E0
                    global base_curr = copy(bij2[j:j+length(w_transformed)-1])
                end
                base_test = base_test[1:end-length(w_transformed)]
            end
            append!(bases, base_curr)
            append!(base_test, base_curr)
            push!(E_list, E_low)
            push!(gaussians, i)
        end
        
        println("Best convergent numerical value: ", E_list[end])
        
    else
        error("Invalid method provided!")
    end
    
    println("Best convergent numerical value: ", E_list[end])
    p = plot(gaussians, E_list, marker=:circle, label="Numerical result",linewidth=2)
    title!(p, "S-wave convergence")
    xlabel!(p, "Number of Gaussians")
    ylabel!(p, "Energy [Hartree]")

    if plot_result==true
        display(p)
    end 
    return p, E_list[end], bases
end
"""
    run_simulation_nuclear(ngauss=2, dim=2, bmax=5)

Run a nuclear simulation and print the final energy.

# Arguments
- `ngauss`: Number of Gaussian functions to use in the simulation (default is 2).
- `dim`: Dimension of the simulation (default is 2).
- `bmax`: Maximum impact parameter (default is 5).

# Outputs
Prints the final energy of the simulation.

# Example
```julia
run_simulation_nuclear(3, 3, 10)
"""
function run_simulation_nuclear(ngauss, dim, bmax, masses, params)
    return E_list, gaussians, eigenvectors, coords, masses = OptimizeGlobalParameters(ngauss, dim, bmax, masses, params)
end

