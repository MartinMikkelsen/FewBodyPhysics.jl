using Plots
export corput, halton, run_simulation, run_simulation_nuclear, plot_convergence

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
function run_simulation(num_gauss=25::Int, method=:quasirandom::Symbol, plot_result=true::Bool)
    b1 = 10
    E_list = []
    gaussians = []
    E_theory = []
    bij = []
    E_S = -0.527 #Ground state energy of hydrogen anion in Hartree
    E_low = Inf
    global bases = []
    global base_test = []
    if method == :quasirandom
        println("---------QUASI-RANDOM METHOD---------")

        for i in 1:num_gauss
            hal = halton(i, 15 * length(w_trans))
            bij = -log.(hal) .* b1
            for j in 1:length(w_trans):length(hal)
                append!(base_test, bij[j:j+length(w_trans)-1])
                E0 = S_energy(base_test, K_trans, w_trans)
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
            push!(E_theory, E_S)
        end

    elseif method == :quasirandomrefined
        println("---------QUASI-RANDOM METHOD W. REFINEMENT---------")

        for i in 1:num_gauss
            hal = halton(i, 15 * length(w_trans))
            bij = -log.(hal) .* b1
            for j in 1:length(w_trans):length(hal)
                append!(base_test, bij[j:j+length(w_trans)-1])
                E0 = S_energy(base_test, K_trans, w_trans)
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
            push!(E_theory, E_S)
        end

        bases_ref = copy(bases)
        E_ref = E_list[end]
        E_list_ref = []
        for i in 1:length(bases_ref)-length(w_trans)
            rand_ref = rand(200 * length(w_trans))
            bij_ref = .-log.(rand_ref) * b1
            for j in 1:length(w_trans):length(rand_ref)
                bases_ref[i:i+length(w_trans)-1] = bij_ref[j:j+length(w_trans)-1]
                E_test = S_energy(bases_ref, K_trans, w_trans)
                if E_test < E_ref
                    E_ref = E_test
                    bases[i:i+length(w_trans)-1] = bij_ref[j:j+length(w_trans)-1]
                end
            end
            bases_ref = copy(bases)
            push!(E_list_ref, E_ref)
            println("E_ref: ", E_ref)
        end

        println("Energy after refinement: ", E_ref)
        println("Difference in energy from before refinement: ", abs(E_ref - E_list[end]))
        println("Difference from target value: ", abs(E_ref - E_S))

    elseif method == :psudorandom
        println("---------PSEUDO-RANDOM METHOD (RANDOM GUESSING)---------")
        for i in 1:num_gauss
            rnd = rand(400 * length(w_trans))
            bij2 = -log.(rnd) * b1
            for j in 1:length(w_trans):length(rnd)
                append!(base_test, bij2[j:j+length(w_trans)-1])
                E0 = S_energy(base_test, K_trans, w_trans)
                if E0 <= E_low
                    E_low = E0
                    global base_curr = copy(bij2[j:j+length(w_trans)-1])
                end
                base_test = base_test[1:end-length(w_trans)]
            end
            append!(bases, base_curr)
            append!(base_test, base_curr)
            push!(E_list, E_low)
            push!(gaussians, i)
            push!(E_theory, E_S)
        end
        
        println("Best convergent numerical value: ", E_list[end])
        println("Theoretical value: ", E_S)
        println("Difference: ", abs(E_list[end] - E_S))
        
    else
        error("Invalid method provided!")
    end
    
    println("Best convergent numerical value: ", E_list[end])
    println("Theoretical value: ", E_S)
    println("Difference: ", abs(E_list[end] - E_S))
    p = plot(gaussians, E_list, marker=:circle, label="Numerical result")
    plot!(p, gaussians, E_theory, linestyle=:dash, label="Theoretical value")
    title!(p, "S-wave convergence of Positron and two Electron System")
    xlabel!(p, "Number of Gaussians")
    ylabel!(p, "Energy [Hartree]")

    if plot_result==true
        display(p)
    end 
    return p
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
function run_simulation_nuclear(ngauss=2, dim=2, bmax=5)
    return E_list, gaussians, eigenvectors, coords, masses = OptimizeGlobalParameters(ngauss, dim, bmax, [mbare, m_π], [b, S])
end

"""
        plot_convergence(plot_params, title, xlabel, ylabel, legendtxt)

Create a scatter plot of the energy list (`"E_list"`) against the number of Gaussians (`"gaussians"`)
from the `plot_params` dictionary. The plot is labeled according to `title`, `xlabel`, `ylabel`, 
and `legendtxt`.

# Arguments
- `plot_params::Dict`: A dictionary containing the data to plot. It should have keys `"gaussians"` 
    and `"E_list"` corresponding to the x and y values of the plot, respectively.
- `title::String`: The title of the plot.
- `xlabel::String`: The label for the x-axis.
- `ylabel::String`: The label for the y-axis.
- `legendtxt::String`: The label for the legend.
"""
function plot_convergence(plot_params, title, xlabel, ylabel, legendtxt)
    plot(plot_params["gaussians"], plot_params["E_list"], marker = :dot, label=legendtxt)
    title!(title)
    xlabel!(xlabel)
    ylabel!(ylabel)
    display(plot())
end