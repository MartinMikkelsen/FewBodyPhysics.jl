using FewBodyPhysics
using Test

@testset "FewBodyPhysics.jl" begin
    # Test Ω function
    @testset "Ω function tests" begin
        masses_1 = [1.0, 2.0, 3.0]
        J, U = Ω(masses_1)
        @test size(J) == (2, 2)
        @test size(U) == (2, 2)
        @test J ≈ [0.0 -0.25; 0.75 0.5] atol=1e-10
        @test U ≈ [1.3333333333 0.0; -0.6666666667 2.0] atol=1e-10

        masses_2 = [1.0]
        J, U = Ω(masses_2)
        @test size(J) == (0, 0)
        @test size(U) == (0, 0)
    end

    # Test A_generate function
    @testset "A_generate function tests" begin
        bij = [1.0, 2.0, 3.0]
        w_list = [rand(3), rand(3), rand(3)]
        A = A_generate(bij, w_list)
    end

    # Test transform_list function
    @testset "transform_list function tests" begin
        α = [1.0, 2.0, 3.0]
        transformed = transform_list(α)
        @test length(transformed) == 3
        @test transformed[1] == [1.0]
        @test transformed[2] == [2.0]
        @test transformed[3] == [3.0]
    end

    # Test shift function
    @testset "shift function tests" begin
        a = [1.0 2.0; 3.0 4.0]
        b = [5.0 6.0; 7.0 8.0]
        mat = [1.0 0.0; 0.0 2.0]
        @test shift(a, b) == 70.0
        @test shift(a, b, mat) == 150.0
    end

    # Test w_gen function
    @testset "w_gen function tests" begin
        dim = 3
        i = 1
        j = 2
        w = w_gen(dim, i, j)
        @test length(w) == dim
        @test w == [1, -1, 0]

        dim = 1
        i = 1
        j = 1
        w = w_gen(dim, i, j)
        @test length(w) == dim
        @test w == [1]
    end

    # Test transform_coordinates and transform_back functions
    @testset "transform_coordinates and transform_back function tests" begin
        masses = [1.0, 2.0, 3.0]
        Ω_matrix, _ = Ω(masses)

        r = [1.0, 2.0, 3.0]
        transformed_r = transform_coordinates(Ω_matrix, r)
        @test length(transformed_r) == 2

        x = [1.0 2.0; 3.0 4.0]
        transformed_x = transform_back(Ω_matrix, x)
        @test size(transformed_x) == (3, 2)
    end
end
