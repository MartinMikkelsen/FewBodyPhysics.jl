var documenterSearchIndex = {"docs":
[{"location":"resources/#Resources","page":"Resources","title":"Resources","text":"","category":"section"},{"location":"resources/","page":"Resources","title":"Resources","text":"Stochastic Variational Approach to Quantum-Mechanical Few-Body Problems, Yasuyuki Suzuki , Kálmán Varga\nThreshold photoproduction of neutral pions off protons in nuclear model with explicit mesons\nCorrelated Gaussians and low-discrepancy sequences\nQuasi-One-Dimensional Few-Body Systems with Correlated Gaussians\nCorrelated Gaussian method in Quantum Mechanics","category":"page"},{"location":"theory/#Theory","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"theory/#Basis-expansion-of-the-Schrödinger-equation","page":"Theory","title":"Basis expansion of the Schrödinger equation","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"We are going to solve the Schrödinger equation","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"hatHpsirangle = epsilonpsirangle","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where hatH is the Hamiltonian of a quantum few-body system, psirangle and epsilon are the eigenfunction and the eigenvalue to be found.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"We shall expand the wave-function psirangle in terms of a set of basis functions irangle for i = 1 ldots n,","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"psirangle = sum_i=1^n c_i irangle","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Inserting the expansion into the Schrödinger equation and multiplying from the left with langle k for 1 leq k leq n gives","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"sum_i=1^n langle khatHirangle c_i = epsilon sum_i=1^n langle kirangle c_i","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Or, in the matrix notation","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"Hc = epsilon Nc","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where H and N are correspondingly the Hamiltonian and the overlap matrices with the matrix elements","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"H_ki = langle khatHirangle quad N_ki","category":"page"},{"location":"theory/#Gaussians-as-basis-functions","page":"Theory","title":"Gaussians as basis functions","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"We shall use the so-called Correlated Gaussians (or Explicitly Correlated Gaussians) as the basis functions. For a system of N particles with coordinates vecr_i, i = 1 ldots N, the Correlated Gaussian is defined as","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"g(vecr_1 ldots vecr_N) = exp left( - sum_ij=1^N A_ijvecr_i cdot vecr_j - sum_i=1^N vecs_i cdot vecr_i right)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where vecr_i cdot vecr_j denotes the dot-product of the two vectors; and where A, a symmetric positive-defined matrix, and vecs_i, i=1ldotsN, the shift-vectors, are (cleverly chosen) parameters of the Gaussian.","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"In matrix notation,","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"g(vecr) = exp left( -vecr^T A vecr + vecs^T vecr right)","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where vecr is the column of the coordinates vecr_i and vecs is the column of the shift-vectors vecs_i,","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"vecr =\nbeginpmatrix\nvecr_1 \nvdots \nvecr_N\nendpmatrix quad\nvecs =\nbeginpmatrix\nvecs_1 \nvdots \nvecs_N\nendpmatrix","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"and","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"vecr^T A vecr + vecs^T vecr = sum_ij vecr_i cdot A_ijvecr_j + sum_i vecs_i cdot vecr_i","category":"page"},{"location":"#FewBodyPhysics.jl","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"","category":"section"},{"location":"#Installation","page":"FewBodyPhysics.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"Get the latest stable release with Julia's package manager:","category":"page"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"julia ] add FewBodyPhysics","category":"page"},{"location":"#Quick-example","page":"FewBodyPhysics.jl","title":"Quick example","text":"","category":"section"},{"location":"","page":"FewBodyPhysics.jl","title":"FewBodyPhysics.jl","text":"using Plots, FewBodyPhysics\n\nw_list = [ [1, -1, 0], [1, 0, -1], [0, 1, -1] ]\nmasses = [1, 1, 1]\nK = [0 0 0; 0 1/2 0; 0 0 1/2]\nJ, U = Ω(masses)\nK_transformed = J * K * J'\nw_transformed = [U' * w_list[i] for i in 1:length(w_list)];\n\np = run_simulation(50, :quasirandom, w_transformed, K_transformed)","category":"page"},{"location":"API/#API","page":"API","title":"API","text":"","category":"section"},{"location":"API/#Coordinates","page":"API","title":"Coordinates","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Ω\nA_generate\ntransform_list\nshift\nw_gen\ntransform_coordinates\ntransform_back","category":"page"},{"location":"API/#FewBodyPhysics.Ω","page":"API","title":"FewBodyPhysics.Ω","text":"Ω(masses::Array)\n\nCalculate the Jacobi transformation matrix J and its inverse U for a system of particles with specified masses.\n\nArguments\n\nmasses::Array: A vector of masses for the particles.\n\nReturns\n\nJ::Matrix: The Jacobi transformation matrix.\nU::Matrix: The inverse of the Jacobi transformation matrix.\n\nNotes\n\nFor systems with more than one particle, the returned matrices exclude the last row/column for proper dimensionality in transformations.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.A_generate","page":"API","title":"FewBodyPhysics.A_generate","text":"A_generate(bij::Array, w_list::Array)\n\nGenerate a matrix A for Gaussian basis functions given width parameters bij and weight vectors w_list.\n\nArguments\n\nbij::Array: A vector of width parameters for the Gaussian basis functions.\nw_list::Array: A list of weight vectors.\n\nReturns\n\nMatrix: The sum of weighted outer products of w_list, scaled by bij.\n\nNotes\n\nThis function is used to construct basis elements for the expansion of few-body wavefunctions.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.transform_list","page":"API","title":"FewBodyPhysics.transform_list","text":"transform_list(α::Array)\n\nTransform a list of scalar values α into a list of 1x1 matrices.\n\nArguments\n\nα::Array: A list of scalar values.\n\nReturns\n\nArray: A list of 1x1 matrices where each matrix contains one of the scalar values from α.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.shift","page":"API","title":"FewBodyPhysics.shift","text":"shift(a::Array, b::Array, mat::Matrix=I)\n\nCalculate the weighted sum of the element-wise product of vectors a and b using matrix mat.\n\nArguments\n\na::Array: A vector or matrix.\nb::Array: A vector or matrix of the same size as a.\nmat::Matrix: An optional matrix to weight the product (default is the identity matrix).\n\nReturns\n\nFloat64: The weighted sum of products.\n\nNotes\n\na and b are typically shift vectors in the configuration space of a few-body system.\nIf mat is provided, its dimensions must match the number of elements in a and b.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.w_gen","page":"API","title":"FewBodyPhysics.w_gen","text":"w_gen(dim::Int, i::Int, j::Int)\n\nGenerate a weight vector for the i-th and j-th coordinates in a space of dimension dim.\n\nArguments\n\ndim::Int: The dimension of the space.\ni::Int: The index for the positive element in the weight vector.\nj::Int: The index for the negative element in the weight vector.\n\nReturns\n\nVector{Int}: A vector with 1 at the i-th position, -1 at the j-th position, and 0 elsewhere.\n\nNotes\n\nThis function is useful for generating basis vectors in few-body coordinate transformations.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.transform_coordinates","page":"API","title":"FewBodyPhysics.transform_coordinates","text":"transform_coordinates(Ω::Matrix{Float64}, r::Vector{Float64})\n\nTransform the coordinates r of a system using the Jacobi matrix Ω.\n\nArguments\n\nΩ::Matrix{Float64}: The Jacobi transformation matrix.\nr::Vector{Float64}: The coordinates to be transformed.\n\nReturns\n\nVector{Float64}: The transformed coordinates.\n\nNotes\n\nThis function applies the inverse of Jacobi matrix J to the coordinate vector r.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.transform_back","page":"API","title":"FewBodyPhysics.transform_back","text":"transform_back(Ω::Matrix{Float64}, x::Matrix{Float64})\n\nTransform the coordinates x back to the original system using the inverse of the Jacobi matrix Ω.\n\nArguments\n\nΩ::Matrix{Float64}: The Jacobi transformation matrix.\nx::Matrix{Float64}: The coordinates to be transformed back.\n\nReturns\n\nMatrix{Float64}: The coordinates transformed back to the original system.\n\nNotes\n\nThis function applies the inverse of matrix U to the coordinate matrix x.\n\n\n\n\n\n","category":"function"},{"location":"API/#Matrix-elements","page":"API","title":"Matrix elements","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"S_elements\nS_wave\nS_energy\nP_elements\npion_nucleon\nComputeEigenSystem\nGetMinimumEnergy\nOptimizeGlobalParameters","category":"page"},{"location":"API/#FewBodyPhysics.S_elements","page":"API","title":"FewBodyPhysics.S_elements","text":"S_elements(A, B, K, w=nothing)\n\nCalculate matrix elements for overlap, kinetic energy, and optionally the Coulomb term.\n\nArguments\n\nA::Matrix: Matrix representing the width of Gaussian basis functions for state i.\nB::Matrix: Matrix representing the width of Gaussian basis functions for state j.\nK::Matrix: Kinetic energy matrix.\nw::Vector (optional): Weight vectors for the particles involved.\n\nReturns\n\nM0::Float64: The overlap matrix element between the two states.\ntra::Float64: The trace used in the kinetic energy calculation.\nCoulomb_term::Float64 (optional): The Coulomb interaction term, if weight vectors w are provided.\n\nNotes\n\nThe Coulomb term is calculated only if the weight vectors w are specified.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.S_wave","page":"API","title":"FewBodyPhysics.S_wave","text":"S_wave(α, K, w=nothing)\n\nCalculate the wavefunction overlap, kinetic energy, and optionally Coulomb interaction matrices for a given set of basis functions.\n\nArguments\n\nα::Vector: A list of scalar width parameters for the Gaussian basis functions.\nK::Matrix: Kinetic energy matrix.\nw::Vector (optional): Weight vectors for the particles involved.\n\nReturns\n\noverlap::Matrix: The overlap matrix for the basis functions.\nkinetic::Matrix: The kinetic energy matrix for the basis functions.\nCoulomb::Matrix (optional): The Coulomb interaction matrix, if weight vectors w are specified.\n\nNotes\n\nThe Coulomb matrix is computed only if the weight vectors w are specified.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.S_energy","page":"API","title":"FewBodyPhysics.S_energy","text":"S_energy(bij, K, w)\n\nCompute the ground state energy of the system using the basis functions specified by the width parameters bij.\n\nArguments\n\nbij::Vector: A list of width parameters for the Gaussian basis functions.\nK::Matrix: Kinetic energy matrix.\nw::Vector: Weight vectors for the particles involved.\n\nReturns\n\nE0::Float64: The lowest eigenvalue computed from the Hamiltonian, considered as the ground state energy of the system.\n\nNotes\n\nThis function constructs the Hamiltonian from the overlap, kinetic, and Coulomb matrices and solves for its eigenvalues.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.P_elements","page":"API","title":"FewBodyPhysics.P_elements","text":"P_elements(a, b, A, B, K, w=nothing)\n\nCalculate the perturbation matrix elements given two basis states represented by vectors a and b, and their respective width matrices A and B.\n\nArguments\n\na::Vector: The coefficient vector for basis state i.\nb::Vector: The coefficient vector for basis state j.\nA::Matrix: Matrix representing the width of Gaussian basis functions for state i.\nB::Matrix: Matrix representing the width of Gaussian basis functions for state j.\nK::Matrix: Kinetic energy matrix.\nw::Vector (optional): Weight vectors for the particles involved.\n\nReturns\n\nM1::Float64: The overlap perturbation term.\nkinetic::Float64: The kinetic energy perturbation term.\nCoulomb_term::Float64 (optional): The Coulomb interaction perturbation term, if weight vectors w are provided.\n\nNotes\n\nThe Coulomb interaction perturbation term is calculated only if the weight vectors w are specified.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.pion_nucleon","page":"API","title":"FewBodyPhysics.pion_nucleon","text":"pion_nucleon(alphas, masses, params)\n\nCalculate the overlap and kinetic matrices for a pion-nucleon system.\n\nArguments\n\nalphas: A vector of alpha values, which are parameters related to the Gaussian basis functions.\nmasses: A 2-element vector containing the masses of the nucleon and the pion, respectively.\nparams: A 2-element vector containing the parameters b and S.\n\nReturns\n\noverlap: A matrix representing the overlap between the basis functions.\nkinetic: A matrix representing the kinetic energy elements.\n\nDescription\n\nThe function calculates the overlap and kinetic matrices for a pion-nucleon system using Gaussian basis functions. The overlap matrix elements are calculated as integrals of the product of two basis functions, and the kinetic matrix elements are calculated as integrals of the product of the derivatives of two basis functions.\n\nThe function uses the reduced mass of the pion-nucleon system, the parameter b related to the width of the Gaussian functions, and the parameter S related to the amplitude of the Gaussian functions.\n\nThe matrices are symmetric, and the diagonal elements of the overlap matrix are 1. The off-diagonal elements are calculated using the alpha parameters and the b and S parameters.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.ComputeEigenSystem","page":"API","title":"FewBodyPhysics.ComputeEigenSystem","text":"ComputeEigenSystem(bs, masses, params)\n\nCalculate the eigenvalues and eigenvectors of a system defined by parameters bs, masses, and params.\n\nArguments\n\nbs: Array of parameter values used in the computation.\nmasses: Array of masses, representing physical properties of the system.\nparams: Additional parameters required for the calculation.\n\nReturns\n\nTuple of eigenvalues (E) and eigenvectors (c).\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.GetMinimumEnergy","page":"API","title":"FewBodyPhysics.GetMinimumEnergy","text":"GetMinimumEnergy(bs, masses, params)\n\nCompute the minimum energy of a system characterized by bs, masses, and params.\n\nArguments\n\nbs: Array of parameter values used in the computation.\nmasses: Array of masses, representing physical properties of the system.\nparams: Additional parameters required for the calculation.\n\nReturns\n\nMinimum energy value of the system.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.OptimizeGlobalParameters","page":"API","title":"FewBodyPhysics.OptimizeGlobalParameters","text":"OptimizeGlobalParameters(ngauss, dim, bmax, masses, params)\n\nPerform global optimization over a given parameter space to find optimal parameters for a physical system.\n\nArguments\n\nngauss: Number of Gaussian functions used in the optimization.\ndim: Dimensionality of the parameter space.\nbmax: Maximum value of the parameter b.\nmasses: Array of masses, representing physical properties of the system.\nparams: Additional parameters used in the optimization.\n\nReturns\n\nE_list: List of optimized energies.\ngaussians: List of Gaussian functions used.\ncoords: Optimized coordinates in the parameter space.\neigenvectors: Eigenvectors corresponding to the optimized coordinates.\nmasses: Updated masses array.\n\n\n\n\n\n","category":"function"},{"location":"API/#Sampling","page":"API","title":"Sampling","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"corput\nhalton\nrun_simulation\nrun_simulation_nuclear","category":"page"},{"location":"API/#FewBodyPhysics.corput","page":"API","title":"FewBodyPhysics.corput","text":"Generate the nth element of the van der Corput sequence in base b.\n\nArguments\n\nn: The nth element of the sequence to be generated.\nb: The base for the van der Corput sequence. Default is 3.\n\nReturns\n\nq: The nth element of the van der Corput sequence in base b.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.halton","page":"API","title":"FewBodyPhysics.halton","text":"Generate the nth d-dimensional point in the Halton sequence.\n\nArguments\n\nn: The nth element of the sequence to be generated.\nd: The dimension of the space.\n\nReturns\n\nAn array containing the nth d-dimensional point in the Halton sequence.\n\nErrors\n\nThrows an assertion error if d exceeds the number of basis elements.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.run_simulation","page":"API","title":"FewBodyPhysics.run_simulation","text":"Run the simulation for a quantum system using quasi-random or pseudo-random methods to determine the S-wave convergence.\n\nArguments\n\nnum_gauss::Int: The number of Gaussians to use in the simulation. Default is 15.\nmethod::Symbol: The method to use for the simulation. Can be :quasirandom, :quasirandomrefined, or :psudorandom. Default is :quasirandom.\nplot_result::Bool: Whether to plot the results. Default is true.\n\nReturns\n\np: The plot object if plot_result is true.\n\nNotes\n\nThe function prints various convergence information and, if plot_result is true, displays a plot of the numerical result against the theoretical value.\n\n\n\n\n\n","category":"function"},{"location":"API/#FewBodyPhysics.run_simulation_nuclear","page":"API","title":"FewBodyPhysics.run_simulation_nuclear","text":"run_simulation_nuclear(ngauss=2, dim=2, bmax=5)\n\nRun a nuclear simulation and print the final energy.\n\nArguments\n\nngauss: Number of Gaussian functions to use in the simulation (default is 2).\ndim: Dimension of the simulation (default is 2).\nbmax: Maximum impact parameter (default is 5).\n\nOutputs\n\nPrints the final energy of the simulation.\n\nExample\n\n```julia runsimulationnuclear(3, 3, 10)\n\n\n\n\n\n","category":"function"}]
}
