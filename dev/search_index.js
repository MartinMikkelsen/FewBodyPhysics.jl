var documenterSearchIndex = {"docs":
[{"location":"resources/#Resources","page":"Resources","title":"Resources","text":"","category":"section"},{"location":"resources/","page":"Resources","title":"Resources","text":"Stochastic Variational Approach to Quantum-Mechanical Few-Body Problems, Yasuyuki Suzuki , Kálmán Varga\nThreshold photoproduction of neutral pions off protons in nuclear model with explicit mesons\nCorrelated Gaussians and low-discrepancy sequences\nQuasi-One-Dimensional Few-Body Systems with Correlated Gaussians\nCorrelated Gaussian method in Quantum Mechanics","category":"page"},{"location":"theory/#Theory","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"theory/#Basis-expansion-of-the-Schrödinger-equation","page":"Theory","title":"Basis expansion of the Schrödinger equation","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"We are going to solve the Schrödinger equation","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"hatHpsirangle = epsilonpsirangle","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where hatH is the Hamiltonian of a quantum few-body system, psirangle and epsilon are the eigenfunction and the eigenvalue to be found.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"We shall expand the wave-function psirangle in terms of a set of basis functions irangle for i = 1 ldots n,","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"psirangle = sum_i=1^n c_i irangle","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Inserting the expansion into the Schrödinger equation and multiplying from the left with langle k for 1 leq k leq n gives","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"sum_i=1^n langle khatHirangle c_i = epsilon sum_i=1^n langle kirangle c_i","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Or, in the matrix notation","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Hc = epsilon Nc","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where H and N are correspondingly the Hamiltonian and the overlap matrices with the matrix elements","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"H_ki = langle khatHirangle quad N_ki","category":"page"},{"location":"theory/#Gaussians-as-basis-functions","page":"Theory","title":"Gaussians as basis functions","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"We shall use the so-called Correlated Gaussians (or Explicitly Correlated Gaussians) as the basis functions. For a system of N particles with coordinates vecr_i, i = 1 ldots N, the Correlated Gaussian is defined as","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"g(vecr_1 ldots vecr_N) = exp left( - sum_ij=1^N A_ijvecr_i cdot vecr_j - sum_i=1^N vecs_i cdot vecr_i right)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where vecr_i cdot vecr_j denotes the dot-product of the two vectors; and where A, a symmetric positive-defined matrix, and vecs_i, i=1ldotsN, the shift-vectors, are (cleverly chosen) parameters of the Gaussian.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"In matrix notation,","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"g(vecr) = exp left( -vecr^T A vecr + vecs^T vecr right)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where vecr is the column of the coordinates vecr_i and vecs is the column of the shift-vectors vecs_i,","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"vecr =\nbeginpmatrix\nvecr_1 \nvdots \nvecr_N\nendpmatrix quad\nvecs =\nbeginpmatrix\nvecs_1 \nvdots \nvecs_N\nendpmatrix","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"and","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"vecr^T A vecr + vecs^T vecr = sum_ij vecr_i cdot A_ijvecr_j + sum_i vecs_i cdot vecr_i","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Suppose you want to calculate the ground state energy of the hydrogen anion in the rest-frame of the proton. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Plots, FewBodyPhysics\n\nmasses = [1e15, 1.0, 1.0]\nw_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]\n\nK = [0 0 0; 0 1/2 0; 0 0 1/2]\nJ, U = jacobi_transform(masses)\nK_transformed = J * K * J'\nw_transformed = [U' * w for w in w_list]\nTheortical_value = -0.527751016523\n\np, Energy, bases = run_simulation(50, :quasirandom, w_transformed, K_transformed)\nplot(p)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"@show Theortical_value = -0.527751016523 - Energy","category":"page"},{"location":"#FewBodyPhysics.jl","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"","category":"section"},{"location":"#Installation","page":"FewBodyPhysics.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"Get the latest stable release with Julia's package manager:","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"julia ] add FewBodyPhysics","category":"page"},{"location":"#Example","page":"FewBodyPhysics.jl","title":"Example","text":"","category":"section"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"We consider a system of a positron and two electrons. The energy of this system has been very accurately calculated by various approaches and it has been found to be -0.262005 in atomic units (a.u.). We calculate the ground-state energy of this systems using correlated Gaussian bases constructed stochastically with pseudorandom and quasirandom sequences. The Hamiltonian of the system is given by","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"H = - sum_i=1^3 frac12m_ifracpartial^2partial boldsymbolr_i^2 + sum_ij=1^3 fracq_i q_jboldsymbolr_i-boldsymbolr_j","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"The masses of the three constituents are m_i = {1, 1, 1} and the charges q_i = {+1, −1, −1}. We can estimate the ground state of this Coulombic three-body system using 50 Gaussians","category":"page"},{"location":"#Quasirandom","page":"FewBodyPhysics.jl","title":"Quasirandom","text":"","category":"section"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"using Plots, FewBodyPhysics\n\nw_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]\nmasses = [1.0,1.0,1.0]\nK = [1/2 0 0; 0 1/2 0; 0 0 1/2]\nJ, U = jacobi_transform(masses)\nK_transformed = J * K * J'\nw_transformed = [U' * w for w in w_list]\n\nTheortical_value = -0.2620050702328\n\np, Energy = run_simulation(50,:quasirandom, w_transformed, K_transformed)\nplot(p)","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"With a difference in energy","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"@show Energy-Theortical_value","category":"page"},{"location":"#Psudorandom","page":"FewBodyPhysics.jl","title":"Psudorandom","text":"","category":"section"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"And similarly for a psudorandom","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"using Plots, FewBodyPhysics\n\nw_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]\nmasses = [1.0,1.0,1.0]\nK = [1/2 0 0; 0 1/2 0; 0 0 1/2]\nJ, U = jacobi_transform(masses)\nK_transformed = J * K * J'\nw_transformed = [U' * w for w in w_list]\n\nTheortical_value = -0.2620050702328\n\np, Energy = run_simulation(50,:psudorandom, w_transformed, K_transformed)\nplot(p)","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"With a difference in energy","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"@show Theortical_value - Energy","category":"page"},{"location":"#Custom-system","page":"FewBodyPhysics.jl","title":"Custom system","text":"","category":"section"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"One of the strengths of this numerical method is that all the matrix elements are analytical. Consider the following example from Threshold photoproduction of neutral pions off protons in nuclear model with explicit mesons. ","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"using Optim, FewBodyPhysics, Plots\n\nb = 3.9\nS = 41.5\n\nparams = [b, S]\nmasses = [(m_p+m_n)/2, m_pi]\n\nenergies, Gaussians, eigenvectors, coordinates, masses = run_simulation_nuclear(5,2,5,masses,params)\n\nrmax = 5 * b\nrmin = 0.01 * b\nstart = log(rmin)\nstop = log(rmax)\ngrid_points = range(rmin,rmax,3000)\n\nΦ = zeros(length(grid_points), length(coordinates))\n\nfor i in eachindex(coordinates)\n    local ϕ = zeros(length(grid_points))\n    ϕ_sum = zeros(length(grid_points))\n    rs = coordinates[i]\n    c = eigenvectors[i]\n    A = [1 / (b^2) for b in rs]\n    for j in 2:min(length(c), length(A)) \n        ϕ_sum .+= c[j] .* exp.(-(A[j-1]) .* grid_points.^2)\n    end\n    ϕ .+= ϕ_sum \n    Φ[:, i] = ϕ ./ c[1]\nend\n\nr = range(rmin,rmax, length=3000)\n\nplot(r, Φ[:,1], title=\"Φ(r)\", label=\"Φ(r)\",ylabel=\"Φ\",xlabel=\"r\", linewidth=2) #with phase\nΦ_prime = diff(Φ[:,1]) ./ diff(r)\nplot!(r[1:end-1], Φ_prime, label=\"Φ'(r)\", linewidth=2)","category":"page"},{"location":"API/#API","page":"API","title":"API","text":"","category":"section"},{"location":"API/#Coordinates","page":"API","title":"Coordinates","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"ParticleSystem\njacobi_transform\ngenerate_A_matrix\ntransform_list\nshift_vectors\ngenerate_weight_vector\ntransform_coordinates\ninverse_transform_coordinates","category":"page"},{"location":"API/#FewBodyPhysics.ParticleSystem","page":"API","title":"FewBodyPhysics.ParticleSystem","text":"ParticleSystem(masses::Vector{Float64})\n\nA data structure representing a system of particles, storing their masses and associated Jacobi transformation matrices for coordinate transformations.\n\nFields\n\nmasses::Vector{Float64}: A vector containing the masses of the particles.\nJ::Matrix{Float64}: The Jacobi transformation matrix, used to convert particle coordinates into Jacobi coordinates.\nU::Matrix{Float64}: The pseudoinverse of the Jacobi transformation matrix J, used to transform Jacobi coordinates back to the particle coordinate system.\n\nConstructor\n\nParticleSystem(masses::Vector{Float64}): Constructs a new ParticleSystem instance.\nArguments:\nmasses: A vector of particle masses. At least two masses are required.\n\n\n\n\n\n","category":"type"},{"location":"API/#FewBodyPhysics.jacobi_transform","page":"API","title":"FewBodyPhysics.jacobi_transform","text":"jacobi_transform(masses::Vector{Float64})::Tuple{Matrix{Float64}, Matrix{Float64}}\n\nCompute the Jacobi transformation matrices J and U for a system of particles with specified masses.\n\nArguments\n\nmasses::Vector{Float64}: A vector of masses for the particles.\n\nReturns\n\n(J::Matrix{Float64}, U::Matrix{Float64}): The Jacobi transformation matrix and its pseudoinverse.\n\nNotes\n\nThe matrices J and U are used to transform between particle coordinates and Jacobi coordinates.\nThe pseudoinverse U is used instead of the inverse to handle cases where J is not square.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.generate_A_matrix","page":"API","title":"FewBodyPhysics.generate_A_matrix","text":"generate_A_matrix(bij::Vector{Float64}, w_list::Vector{Vector{Float64}})::Matrix{Float64}\n\nGenerate the matrix A for Gaussian basis functions given width parameters bij and weight vectors w_list.\n\nArguments\n\nbij::Vector{Float64}: A vector of width parameters for the Gaussian basis functions.\nw_list::Vector{Vector{Float64}}: A list of weight vectors.\n\nReturns\n\nA::Matrix{Float64}: The sum of weighted outer products of w_list, scaled by bij.\n\nNotes\n\nThis function constructs the A matrix used in the correlated Gaussian method.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.transform_list","page":"API","title":"FewBodyPhysics.transform_list","text":"transform_list(α::Vector{Float64})::Vector{Matrix{Float64}}\n\nTransform a list of scalar values α into a list of 1x1 matrices.\n\nArguments\n\nα::Vector{Float64}: A list of scalar values.\n\nReturns\n\nArray{Matrix{Float64}}: A list of 1x1 matrices where each matrix contains one of the scalar values from α.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.shift_vectors","page":"API","title":"FewBodyPhysics.shift_vectors","text":"shift_vectors(a::Matrix{Float64}, b::Matrix{Float64}, mat::Union{Nothing, Matrix{Float64}}=nothing)::Float64\n\nCalculate the weighted sum of the element-wise product of vectors a and b using matrix mat.\n\nArguments\n\na::Matrix{Float64}: A matrix where each column is a vector a_i.\nb::Matrix{Float64}: A matrix where each column is a vector b_j.\nmat::Union{Nothing, Matrix{Float64}}: An optional matrix to weight the product (default is the identity matrix).\n\nReturns\n\nsum_val::Float64: The weighted sum of products.\n\nNotes\n\nThe matrices a and b should have the same dimensions.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.generate_weight_vector","page":"API","title":"FewBodyPhysics.generate_weight_vector","text":"generate_weight_vector(dim::Int, i::Int, j::Int)::Vector{Int}\n\nGenerate a weight vector for the i-th and j-th coordinates in a space of dimension dim.\n\nArguments\n\ndim::Int: The dimension of the space.\ni::Int: The index for the positive element in the weight vector.\nj::Int: The index for the negative element in the weight vector.\n\nReturns\n\nw::Vector{Int}: A vector with 1 at the i-th position, -1 at the j-th position, and 0 elsewhere.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.transform_coordinates","page":"API","title":"FewBodyPhysics.transform_coordinates","text":"transform_coordinates(J::Matrix{Float64}, r::Vector{Float64})::Vector{Float64}\n\nTransform the coordinates r of a system using the Jacobi matrix J.\n\nArguments\n\nJ::Matrix{Float64}: The Jacobi transformation matrix.\nr::Vector{Float64}: The coordinates to be transformed.\n\nReturns\n\nx::Vector{Float64}: The transformed coordinates.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.inverse_transform_coordinates","page":"API","title":"FewBodyPhysics.inverse_transform_coordinates","text":"inverse_transform_coordinates(U::Matrix{Float64}, x::Vector{Float64})::Vector{Float64}\n\nTransform the coordinates x back to the original system using the inverse of the Jacobi matrix.\n\nArguments\n\nU::Matrix{Float64}: The inverse Jacobi transformation matrix.\nx::Vector{Float64}: The coordinates in Jacobi space.\n\nReturns\n\nr::Vector{Float64}: The coordinates transformed back to the original system.\n\n\n\n\n\n","category":"function"},{"location":"API/#Matrix-elements","page":"API","title":"Matrix elements","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"S_elements\nS_wave\nS_energy\nP_elements\npion_nucleon\nComputeEigenSystem\nGetMinimumEnergy\nOptimizeGlobalParameters","category":"page"},{"location":"API/#FewBodyPhysics.S_elements","page":"API","title":"FewBodyPhysics.S_elements","text":"S_elements(A::Matrix{Float64}, B::Matrix{Float64}, K::Matrix{Float64}, w::Union{Nothing, Vector{Vector{Float64}}}=nothing)\n\nCalculate matrix elements for overlap, kinetic energy, and optionally the Coulomb term.\n\nArguments\n\nA: Width matrix of Gaussian basis functions for state i.\nB: Width matrix of Gaussian basis functions for state j.\nK: Kinetic energy matrix.\nw (optional): List of weight vectors for the particles involved.\n\nReturns\n\nM0: The overlap matrix element between the two states.\ntrace: The trace used in the kinetic energy calculation.\nCoulomb_term (optional): The Coulomb interaction term, if weight vectors w are provided.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.S_wave","page":"API","title":"FewBodyPhysics.S_wave","text":"S_wave(α::Vector{Matrix{Float64}}, K::Matrix{Float64}, w::Union{Nothing, Vector{Vector{Float64}}}=nothing)\n\nCalculate the overlap, kinetic energy, and optionally Coulomb interaction matrices for a given set of basis functions.\n\nArguments\n\nα: A list of width matrices for the Gaussian basis functions.\nK: Kinetic energy matrix.\nw (optional): List of weight vectors for the particles involved.\n\nReturns\n\noverlap: The overlap matrix for the basis functions.\nkinetic: The kinetic energy matrix for the basis functions.\nCoulomb: The Coulomb interaction matrix, if weight vectors w are specified.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.S_energy","page":"API","title":"FewBodyPhysics.S_energy","text":"S_energy(bij::Vector{Float64}, K::Matrix{Float64}, w::Vector{Vector{Float64}})\n\nCompute the ground state energy of the system using the basis functions specified by the width parameters bij.\n\nArguments\n\nbij: A vector of width parameters for the Gaussian basis functions.\nK: Kinetic energy matrix.\nw: List of weight vectors for the particles involved.\n\nReturns\n\nE0: The lowest eigenvalue computed from the Hamiltonian.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.P_elements","page":"API","title":"FewBodyPhysics.P_elements","text":"P_elements(a::Vector{Float64}, b::Vector{Float64}, A::PositiveDefiniteSymmetricMatrix, B::PositiveDefiniteSymmetricMatrix, K::Matrix{Float64}, w::Union{Nothing, Vector{Vector{Float64}}}=nothing)\n\nCalculate the perturbation matrix elements given two basis states represented by vectors a and b, and their respective width matrices A and B.\n\nArguments\n\na: The coefficient vector for basis state i.\nb: The coefficient vector for basis state j.\nA: Width matrix for state i.\nB: Width matrix for state j.\nK: Kinetic energy matrix.\nw (optional): List of weight vectors for the particles involved.\n\nReturns\n\nM1: The overlap perturbation term.\nkinetic: The kinetic energy perturbation term.\nCoulomb_term (optional): The Coulomb interaction perturbation term, if weight vectors w are provided.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.pion_nucleon","page":"API","title":"FewBodyPhysics.pion_nucleon","text":"pion_nucleon(alphas, masses, params)\n\nCalculate the overlap and kinetic matrices for a pion-nucleon system.\n\nArguments\n\nalphas: A vector of alpha values, which are parameters related to the Gaussian basis functions.\nmasses: A 2-element vector containing the masses of the nucleon and the pion, respectively.\nparams: A 2-element vector containing the parameters b and S.\n\nReturns\n\noverlap: A matrix representing the overlap between the basis functions.\nkinetic: A matrix representing the kinetic energy elements.\n\nDescription\n\nThe function calculates the overlap and kinetic matrices for a pion-nucleon system using Gaussian basis functions. The overlap matrix elements are calculated as integrals of the product of two basis functions, and the kinetic matrix elements are calculated as integrals of the product of the derivatives of two basis functions.\n\nThe function uses the reduced mass of the pion-nucleon system, the parameter b related to the width of the Gaussian functions, and the parameter S related to the amplitude of the Gaussian functions.\n\nThe matrices are symmetric, and the diagonal elements of the overlap matrix are 1. The off-diagonal elements are calculated using the alpha parameters and the b and S parameters.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.ComputeEigenSystem","page":"API","title":"FewBodyPhysics.ComputeEigenSystem","text":"ComputeEigenSystem(bs, masses, params)\n\nCalculate the eigenvalues and eigenvectors of a system defined by parameters bs, masses, and params.\n\nArguments\n\nbs: Array of parameter values used in the computation.\nmasses: Array of masses, representing physical properties of the system.\nparams: Additional parameters required for the calculation.\n\nReturns\n\nTuple of eigenvalues (E) and eigenvectors (c).\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.GetMinimumEnergy","page":"API","title":"FewBodyPhysics.GetMinimumEnergy","text":"GetMinimumEnergy(bs, masses, params)\n\nCompute the minimum energy of a system characterized by bs, masses, and params.\n\nArguments\n\nbs: Array of parameter values used in the computation.\nmasses: Array of masses, representing physical properties of the system.\nparams: Additional parameters required for the calculation.\n\nReturns\n\nMinimum energy value of the system.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.OptimizeGlobalParameters","page":"API","title":"FewBodyPhysics.OptimizeGlobalParameters","text":"OptimizeGlobalParameters(ngauss, dim, bmax, masses, params)\n\nPerform global optimization over a given parameter space to find optimal parameters for a physical system.\n\nArguments\n\nngauss: Number of Gaussian functions used in the optimization.\ndim: Dimensionality of the parameter space.\nbmax: Maximum value of the parameter b.\nmasses: Array of masses, representing physical properties of the system.\nparams: Additional parameters used in the optimization.\n\nReturns\n\nE_list: List of optimized energies.\ngaussians: List of Gaussian functions used.\ncoords: Optimized coordinates in the parameter space.\neigenvectors: Eigenvectors corresponding to the optimized coordinates.\nmasses: Updated masses array.\n\n\n\n\n\n","category":"function"},{"location":"API/#Sampling","page":"API","title":"Sampling","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"corput\nhalton\ncalculate_energies\nrun_simulation\nrun_simulation_nuclear","category":"page"},{"location":"API/#FewBodyPhysics.corput","page":"API","title":"FewBodyPhysics.corput","text":"corput(n::Int, b::Int=3) -> Float64\n\nGenerates the nth term of the van der Corput sequence in the given base b, which is often used for quasi-random number generation.\n\nArguments\n\nn::Int: The position of the term in the sequence to calculate.\nb::Int=3: The base for the sequence. Defaults to 3 if not provided.\n\nReturns\n\nFloat64: The nth term in the van der Corput sequence for the given base b.\n\nExample\n\n```julia corput(1, 2)  # Returns 0.5 for base 2\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.halton","page":"API","title":"FewBodyPhysics.halton","text":"halton(n::Int, d::Int) -> Vector{Float64}\n\nGenerates a point in the Halton sequence with d dimensions, used in quasi-random sampling.\n\n# Arguments\n- `n::Int`: The index of the sequence point to generate.  \n- `d::Int`: The dimensionality of the Halton sequence (i.e., the number of bases to use).\n# Returns\nVector{Float64}: A vector of length d representing the nth point in the Halton sequence.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.calculate_energies","page":"API","title":"FewBodyPhysics.calculate_energies","text":"calculate_energies(num_gauss::Int, w_transformed::Vector{Vector{Float64}}, K_transformed::Matrix{Float64}, b1::Float64, method::Symbol) -> Tuple{Vector{Float64}, Vector{Vector{Float64}}, Vector{Int}, Float64}\n\nCalculates and refines energies for a set of Gaussian basis functions using quasi-random or pseudo-random methods.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.run_simulation","page":"API","title":"FewBodyPhysics.run_simulation","text":"run_simulation(num_gauss::Int, method::Symbol, w_transformed::Vector{Vector{Float64}}, K_transformed::Matrix{Float64}, plot_result::Bool=true) -> Tuple{Plots.Plot, Float64, Vector{Vector{Float64}}}\n\nRuns a simulation to calculate energy values for Gaussian basis functions and optionally plots the results.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.run_simulation_nuclear","page":"API","title":"FewBodyPhysics.run_simulation_nuclear","text":"Run a nuclear simulation and print the final energy.\n\n\n\n\n\n","category":"function"},{"location":"API/#Constants","page":"API","title":"Constants","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"μ","category":"page"},{"location":"API/#FewBodyPhysics.μ","page":"API","title":"FewBodyPhysics.μ","text":"μ = (mbare * mpi0) / (mbare + mpi0)\n\nReduced mass of the nucleon-pion system in MeV/c^2.\n\n\n\n\n\n","category":"constant"}]
}