using FewBodyPhysics
using Test
using LinearAlgebra

@testset "FewBodyPhysics.jl" begin
    # Test for Ω function
    @testset "Ω Function Tests" begin
        masses = [1.0, 2.0, 3.0]
        J, U = Ω(masses)
        @test size(J) == (2, 3)
        @test size(U) == (3, 2)
        @test J * U ≈ I(2) atol=1e-10
    end

    # Test for A_generate function
    @testset "A_generate Function Tests" begin
        bij = [1.0, 2.0, 3.0]
        w_list = [rand(3) for _ in 1:3]
        A = A_generate(bij, w_list)
        @test size(A) == (3, 3)
        # test properties of matrix A (e.g., symmetry, trace)
    end

    # Test for transform_list function
    @testset "transform_list Function Tests" begin
        α = [1.0, 2.0, 3.0]
        transformed_list = transform_list(α)
        @test all([isequal(m, [a]) for (m, a) in zip(transformed_list, α)])
    end

    # Test for shift function
    @testset "Shift Function Tests" begin
        a = rand(3)
        b = rand(3)
        mat = rand(3, 3)
        s = shift(a, b, mat)
        @test typeof(s) == Float64
        # More specific tests on the output based on known input
    end

    # Test for w_gen function
    @testset "w_gen Function Tests" begin
        dim = 5
        i = 2
        j = 4
        w = w_gen(dim, i, j)
        @test length(w) == dim
        @test w[i] == 1 && w[j] == -1
        @test sum(w) == 0
    end

    # Test for transform_coordinates function
    @testset "transform_coordinates Function Tests" begin
        masses = [1.0, 2.0, 3.0]
        r = rand(3)
        Ω_matrix = Ω(masses)[1]
        x = transform_coordinates(Ω_matrix, r)
        @test length(x) == length(r) - 1 # as the last dimension is reduced
    end

    # Test for transform_back function
    @testset "transform_back Function Tests" begin
        masses = [1.0, 2.0, 3.0]
        x = rand(3, 2)
        Ω_matrix = Ω(masses)[1]
        r_back = transform_back(Ω_matrix, x)
        @test size(r_back) == size(x)
    end
end
