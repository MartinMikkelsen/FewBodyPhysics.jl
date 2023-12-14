using Plots

b = 3.9
S = 41.5

params = [b, S]
masses = [(m_p+m_n)/2, m_π]

energies, Gaussians, eigenvectors, coordinates, masses = run_simulation_nuclear(5,2,5,masses,params)

rmax = 5 * b
rmin = 0.01 * b
start = log(rmin)
stop = log(rmax)
grid_points = range(rmin,rmax,3000)

Φ = zeros(length(grid_points), length(coordinates))

for i in 1:length(coordinates)
    local ϕ = zeros(length(grid_points))
    ϕ_sum = zeros(length(grid_points))
    rs = coordinates[i]
    c = eigenvectors[i]
    A = [1 / (b^2) for b in rs]
    for j in 2:min(length(c), length(A)) 
        ϕ_sum .+= c[j] .* exp.(-(A[j-1]) .* grid_points.^2)
    end
    ϕ .+= ϕ_sum 
    Φ[:, i] = ϕ ./ c[1]
end

r = range(rmin,rmax, length=3000)

p1 = plot(r, Φ[:,1], title="Φ(r)", label="Φ(r)",ylabel="Φ",xlabel="r", linewidth=2) #with phase
Φ_prime = diff(Φ[:,1]) ./ diff(r)
plot!(r[1:end-1], Φ_prime, label="Φ'(r)", linewidth=2)
p2 = plot(Gaussians, energies, title="Energy = $(round(energies[end]; digits=3))", label="Convergence",ylabel="Energy",xlabel="Number of Gaussians",linewidth=2)
plot(p1, p2, layout = (2, 1))