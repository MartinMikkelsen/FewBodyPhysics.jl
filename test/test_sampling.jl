using Test
using FewBodyPhysics.Sampling

@testset "Sampling Module Tests" begin
    @testset "Van der Corput Sequence" begin
        # Test base 2 (default)
        @test corput(1, 2) ≈ 0.5
        @test corput(2, 2) ≈ 0.25
        @test corput(3, 2) ≈ 0.75
        @test corput(4, 2) ≈ 0.125
        
        # Test with default base parameter
        @test corput(1) ≈ 0.5
        @test corput(2) ≈ 0.25
        
        # Test base 3
        @test corput(1, 3) ≈ 1/3
        @test corput(2, 3) ≈ 2/3
        @test corput(3, 3) ≈ 1/9
        @test corput(4, 3) ≈ 4/9
        
        # Test base 10
        @test corput(1, 10) ≈ 0.1
        @test corput(12, 10) ≈ 0.21
        
        # Test n = 0
        @test corput(0) ≈ 0.0
        
    end
    
    @testset "Halton Sequence" begin
        # Test 1D and 2D sequences
        @test halton(1, 1) ≈ [0.5]
        @test halton(1, 2) ≈ [0.5, 1/3]
        @test halton(2, 2) ≈ [0.25, 2/3]
        
        # Test higher dimensions
        h3 = halton(5, 3)
        @test length(h3) == 3
        @test h3[1] == corput(5, 2)  # First dimension uses base 2
        @test h3[2] == corput(5, 3)  # Second dimension uses base 3
        @test h3[3] == corput(5, 5)  # Third dimension uses base 5
        
        # Test dimension bounds
        @test_throws AssertionError halton(1, 100)
    end
    
    @testset "Bij Generation" begin
        # Test quasirandom generation
        b1 = 2.0
        bij_quasi = generate_bij(:quasirandom, 3, 4, b1)
        @test length(bij_quasi) == 4
        @test all(0 .<= bij_quasi .<= b1)
        @test bij_quasi == halton(3, 4) * b1
        
        # Test invalid method
        @test_throws ErrorException generate_bij(:invalid, 3, 4, b1)
    end
end