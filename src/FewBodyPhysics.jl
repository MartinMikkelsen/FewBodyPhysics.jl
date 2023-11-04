module FewBodyPhysics
using LinearAlgebra, Plots

include("coordinates.jl")
include("matrix_elements.jl")

function run_simulation(num_gauss=15::Int, method=:quasirandom::Symbol)
    b1 = 7
    E_list = []
    gaussians = []
    E_theoS = []
    bij = []
    E_S = -0.527
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

    elseif method == :quasirandomrefined
        println("---------QUASI-RANDOM METHOD W. REFINEMENT---------")

        for i in 1:num_gauss
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

        bases_ref = copy(bases)
        E_ref = E_list[end]
        E_list_ref = []
        for i in 1:length(bases_ref)-length(w_trans)
            rand_ref = rand(200 * length(w_trans))
            bij_ref = .-log.(rand_ref) * b1
            for j in 1:length(w_trans):length(rand_ref)
                bases_ref[i:i+length(w_trans)-1] = bij_ref[j:j+length(w_trans)-1]
                E_test = energyS(bases_ref, K_trans, w_trans)
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
                E0 = energyS(base_test, K_trans, w_trans)
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
            push!(E_theoS, E_S)
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

    plot(gaussians, E_list, marker=:circle, label="Numerical result")
    plot!(gaussians, E_theoS, linestyle=:dash, label="Theoretical value")
    title!("S-wave convergence of Positron and two Electron System")
    xlabel!("Number of Gaussians")
    ylabel!("Energy [Hartree]")


end
run_simulation(15,:quasirandom)
end
