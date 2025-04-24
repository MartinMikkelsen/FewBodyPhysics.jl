using Test
using FewBodyPhysics.Coordinates
using LinearAlgebra
using FewBodyPhysics.Types

@testset "Coordinates Module Tests" begin
    @testset "jacobi_transform" begin
        # Test with simple equal masses
        masses = [1.0, 1.0, 1.0]
        J, U = jacobi_transform(masses)
        
        @test size(J) == (2, 3)
        @test size(U) == (3, 2)
        @test J * U ≈ I(2) atol=1e-10
        
        # Test with different masses
        masses = [1.0, 2.0, 3.0]
        J, U = jacobi_transform(masses)
        
        @test size(J) == (2, 3)
        @test size(U) == (3, 2)
        @test J * U ≈ I(2) atol=1e-10
        
        # Test with exactly two masses
        masses = [1.0, 2.0]
        J, U = jacobi_transform(masses)
        
        @test size(J) == (1, 2)
        @test size(U) == (2, 1)
        @test J * U ≈ [1.0] atol=1e-10
        
        # Test error for fewer than two masses
        @test_throws AssertionError jacobi_transform([1.0])
    end
    
    @testset "ParticleSystem Constructor" begin
        # Test valid construction
        masses = [1.0, 2.0, 3.0]
        ps = ParticleSystem(masses)
        
        @test ps.masses == masses
        @test size(ps.J) == (2, 3)
        @test size(ps.U) == (3, 2)
        @test ps.scale === nothing
        
        # Test with scale
        ps = ParticleSystem(masses, scale=:atomic)
        @test ps.scale === :atomic
        
        # Test error for invalid masses
        @test_throws AssertionError ParticleSystem([1.0])
    end
    
    @testset "default_b0" begin
        @test default_b0(:atomic) == 1.0
        @test default_b0(:molecular) == 3.0
        @test default_b0(:nuclear) == 0.03
        @test default_b0(nothing) == 10.0
        @test_throws ErrorException default_b0(:unknown)
    end
    
    @testset "generate_A_matrix" begin
        bij = [1.0, 2.0]
        w_list = [[1.0, -1.0, 0.0], [0.0, 1.0, -1.0]]
        A = generate_A_matrix(bij, w_list)
        
        @test size(A) == (3, 3)
        @test A[1,1] ≈ 1.0 atol=1e-10
        @test A[2,2] ≈ 1.0 + 0.25 atol=1e-10
        
        # Test error for mismatched lengths
        @test_throws AssertionError generate_A_matrix([1.0], w_list)
        @test_throws AssertionError generate_A_matrix(bij, [[1.0, -1.0], [0.0, 1.0]])
    end
    
    @testset "transform_list" begin
        α = [1.0, 2.0, 3.0]
        result = transform_list(α)
        
        @test length(result) == 3
        @test result[1] == [1.0]
        @test result[2] == [2.0]
        @test result[3] == [3.0]
    end
    
    @testset "shift_vectors" begin
        a = [1.0 2.0; 3.0 4.0]
        b = [5.0 6.0; 7.0 8.0]
        
        # Test with default identity matrix
        result = shift_vectors(a, b)
        expected = dot([1.0, 3.0], [5.0, 7.0]) + dot([2.0, 4.0], [6.0, 8.0])
        @test result ≈ expected atol=1e-10
        
        # Test with custom weighting matrix
        mat = [2.0 1.0; 1.0 3.0]
        result = shift_vectors(a, b, mat)
        expected = 2*dot([1.0, 3.0], [5.0, 7.0]) + dot([1.0, 3.0], [6.0, 8.0]) +
                   dot([2.0, 4.0], [5.0, 7.0]) + 3*dot([2.0, 4.0], [6.0, 8.0])
        @test result ≈ expected atol=1e-10
        
        # Test error for mismatched dimensions
        @test_throws AssertionError shift_vectors(a, [1.0 2.0])
        @test_throws AssertionError shift_vectors(a, b, [1.0 2.0])
    end
    
    @testset "generate_weight_vector" begin
        w = generate_weight_vector(3, 1, 2)
        @test w == [1, -1, 0]
        
        w = generate_weight_vector(4, 2, 4)
        @test w == [0, 1, 0, -1]
        
        @test_throws AssertionError generate_weight_vector(3, 0, 2)
        @test_throws AssertionError generate_weight_vector(3, 1, 4)
    end
    
    @testset "Coordinate transformations" begin
        # Create test data
        masses = [1.0, 2.0, 3.0]
        J, U = jacobi_transform(masses)
        r = [1.0, 2.0, 3.0]
        
        # Test transform_coordinates
        x = transform_coordinates(J, r)
        @test size(x) == (2,)
        
        # Test inverse_transform_coordinates
        r_recovered = inverse_transform_coordinates(U, x)
        @test size(r_recovered) == (3,)
        
        # Test that the round-trip transformation recovers the original coordinates
        @test r_recovered ≈ r atol=1e-10
        
        # Test error handling
        @test_throws AssertionError transform_coordinates(J, [1.0, 2.0])
        @test_throws AssertionError inverse_transform_coordinates(U, [1.0])
    end
end