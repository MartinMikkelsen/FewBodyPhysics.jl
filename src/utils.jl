module Utils

export ψ₀, plot_wavefunction, plot_density

"""
    ψ₀(r::Vector{Float64}, c₀::Vector{Float64}, basis_fns::Vector{<:GaussianBase})

Evaluates the ground state wavefunction ψ₀ at position `r`.
"""
function ψ₀(r::Vector{Float64}, c₀::Vector{Float64}, basis_fns::Vector)
    sum(c₀[i] * exp(-r' * basis_fns[i].A * r) for i in eachindex(basis_fns))
end

"""
    plot_wavefunction(c₀, basis_fns; axis=1, range=(-4.0, 4.0), N=400)

Plots a 1D slice of the wavefunction ψ₀(x, 0) or ψ₀(0, y).
"""
function plot_wavefunction(c₀, basis_fns; axis=1, range=(-4.0, 4.0), N=400)
    xs = range(range[1], range[2], length=N)
    ψ_vals = [ψ₀(axis == 1 ? [x, 0.0] : [0.0, x], c₀, basis_fns) for x in xs]
    plot(xs, real(ψ_vals), xlabel=axis == 1 ? "x" : "y", ylabel="ψ₀",
         lw=2, title="Ground-state Wavefunction", legend=false)
end

"""
    plot_density(c₀, basis_fns; range=(-4.0, 4.0), N=200)

Plots the 2D probability density |ψ₀(x, y)|² in Jacobi space.
"""
function plot_density(c₀, basis_fns; range=(-4.0, 4.0), N=200)
    xs = range(range[1], range[2], length=N)
    ys = range(range[1], range[2], length=N)
    ψ² = [abs2(ψ₀([x, y], c₀, basis_fns)) for y in ys, x in xs]
    heatmap(xs, ys, ψ², xlabel="x", ylabel="y", title="|ψ₀(x, y)|²", aspect_ratio=1)
end

end
