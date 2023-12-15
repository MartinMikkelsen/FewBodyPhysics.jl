using FewBodyPhysics
using Test
using LinearAlgebra

@testset "FewBodyPhysics.jl" begin
    @testset "Ω tests" begin
        J, U = Ω([1, 2, 3])
        @test J isa Matrix
        @test U isa Matrix
        @test size(J) == (2, 3)
        @test size(U) == (3, 2)
    end
    
    @testset "A_generate tests" begin
        A = A_generate([1, 2, 3], [[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        @test A isa Matrix
        @test size(A) == (3, 3)
    end
    
    @testset "transform_list tests" begin
        transformed = transform_list([1, 2, 3])
        @test transformed isa Array
        @test all([isa(x, Matrix) for x in transformed])
        @test all([size(x) == (1, 1) for x in transformed])
    end
    
    @testset "shift tests" begin
        s = shift([1, 2, 3], [4, 5, 6])
        @test s isa Float64
        @test_throws AssertionError shift([1, 2, 3], [4, 5, 6], [1, 2])
    end
    
    @testset "w_gen tests" begin
        w = w_gen(3, 1, 2)
        @test w isa Vector{Int}
        @test size(w) == (3,)
        @test w == [1, -1, 0]
    end
    
    @testset "transform_coordinates tests" begin
        r_transformed = transform_coordinates([1 0; 0 1], [1, 2])
        @test r_transformed isa Vector{Float64}
        @test size(r_transformed) == (2,)
    end
    
    @testset "transform_back tests" begin
        x_transformed = transform_back([1 0; 0 1], [1 2; 3 4])
        @test x_transformed isa Matrix{Float64}
        @test size(x_transformed) == (2, 2)
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
