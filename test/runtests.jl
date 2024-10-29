using FewBodyPhysics
using Test
using LinearAlgebra

@testset "FewBodyPhysics.jl" begin
    @testset "Coordinate Transformations" begin
        masses = [1.0, 2.0, 3.0]
        J, U = JacobiTransform(masses)
        @test size(J) == (2, 3)
        @test size(U) == (3, 2)
    
        bij = [0.5, 0.8]
        w_list = [generate_weight_vector(3, 1, 2), generate_weight_vector(3, 2, 3)]
        A = generate_A_matrix(bij, w_list)
        @test size(A) == (3, 3)
    
        α = [1.0, 2.0, 3.0]
        transformed_list = transform_list(α)
        @test length(transformed_list) == 3
        @test transformed_list[1] == [1.0]
    
        a = rand(3, 2)
        b = rand(3, 2)
        sum_val = shift_vectors(a, b)
        @test isa(sum_val, Float64)
    
        w = generate_weight_vector(3, 1, 2)
        @test w == [1, -1, 0]
    
        r = rand(3)
        x = transform_coordinates(J, r)
        r_back = inverse_transform_coordinates(U, x)
        @test norm(r - r_back) < 1e-10
    end

    @testset "corput tests" begin
        @test corput(1, 2) ≈ 0.5
        @test corput(2, 2) ≈ 0.25
        @test corput(3, 2) ≈ 0.75
    end
    
    @testset "halton tests" begin
        @test halton(1, 2) ≈ [0.5, 0.3333333333333333]
        @test halton(2, 2) ≈ [0.25, 0.6666666666666666]
        @test_throws AssertionError halton(1, 100)
    end
    
    @testset "run_simulation tests" begin
        w_transformed = [1, 2, 3]
        K_transformed = [1, 2, 3]
        p, E_list, bases = run_simulation(15, :quasirandom, w_transformed, K_transformed, false)
        @test E_list isa Float64
        @test bases isa Array
        @test_throws ErrorException run_simulation(15, :invalidmethod, w_transformed, K_transformed, false)
    end
    
    @testset "run_simulation_nuclear tests" begin
        masses = [1, 2, 3]
        params = [1, 2, 3]
        E_list, gaussians, eigenvectors, coords, masses = run_simulation_nuclear(2, 2, 5, masses, params)
        @test E_list isa Array
        @test gaussians isa Array
        @test eigenvectors isa Array
        @test coords isa Array
        @test masses isa Array
    end

    
end
